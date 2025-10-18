
/* attempt to solve this problem using simulated annealing.
 * This is actually kind of a difficult problem to solve.
 * Mostly because of the dependency constraint on previous
 * values. That prevents us from breaking the system into
 * multiple independent equations, and also prevents us
 * from fixing the variable count in our error functions,
 * which makes it difficult to use a lot of existing techniques. */

#define GRAVITY 9.807

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>

struct data_sample {
   double sec;
   double thrust;
};

static void compute_matching_M(struct data_sample * samples, double * M, double Ve, int count) {
   double accum = 0.0; /* previous mass burned */

   for (int i = 0; i < count; ++i) {

      M[i] = (samples[i].thrust + (accum*GRAVITY)) / Ve;

      accum += M[i] * samples[i].sec;

   }
}

/* error */
static double E(struct data_sample * samples, double * M, double Ve, double mB, int count) {

   /* fma is bae... xD */

   double mB_err = -mB;

   for (int i = 0; i < count; ++i)
      mB_err = fma(M[i],samples[i].sec,mB_err);

   mB_err = mB_err * mB_err; /* squared error */

   double err = mB_err;

   /*
   double accum = 0.0; // previous mass burned
   for (int i = 0; i < count; ++i) {

      double f_err = -samples[i].thrust;
      f_err = fma(M[i],Ve,f_err);
      f_err = fma(accum,-GRAVITY,f_err);

      accum += M[i] * samples[i].sec;

      err += f_err * f_err; // squared error
   }
   */

   return err;
}

#define SWAP(x,y)        \
   do {                  \
      double * temp = x; \
      x = y;             \
      y = temp;          \
   } while (0)


int main(int argc, char ** argv) {
   if (argc != 4) {
      printf("usage:\n%s <startN> <endN> <file>\n",argv[0]);
      exit(1);
   }

   int start = strtol(argv[1],NULL,0);
   int end = strtol(argv[2],NULL,0);

   printf("sample range is %d -> %d\n",start,end);
   printf("file name is: %s\n",argv[3]);

   printf("collecting values...\n");

   FILE * fptr = fopen(argv[3],"r");

   if (fptr == NULL) {
      printf("unable to open file %s\n",argv[3]);
      perror("fopen");
      exit(1);
   }

   char buf[1024];
   int epos;

   struct data_sample * data_vec =
                        malloc(sizeof(struct data_sample)*256);
   int data_vec_pos = 0, data_vec_len = 256;

   double prev_sec = 0.0, temp;
   while (fgets(buf,sizeof(buf),fptr)) {

      /* skip any lines where a digit is not the first character */
      if (!isdigit(buf[0])) {
         buf[strlen(buf)-1] = '\0';
         printf("skipping line: \"%s\"\n",buf);
         continue;
      }

      char * pos = buf;
      long collect_time = strtol(pos,&pos,0);
      long raw = strtol(pos,&pos,0);
      double value = strtod(pos,&pos);
      double value_tare = strtod(pos,&pos);

      if (data_vec_pos == data_vec_len) {
         data_vec_len *= 2;
         data_vec =
            realloc(data_vec,sizeof(struct data_sample)*data_vec_len);
      }

      data_vec[data_vec_pos].sec = (double)collect_time / (double)1e6;
      temp = data_vec[data_vec_pos].sec;
      data_vec[data_vec_pos].sec -= prev_sec;
      prev_sec = temp;
      data_vec[data_vec_pos].thrust = value_tare;

      data_vec_pos += 1;
   }

   fclose(fptr);

   if (!data_vec_pos) {
      printf("no entries in file, exiting early\n");
      exit(0);
   }

   if (data_vec_pos <= end) {
      printf("entries <= end -> %d <= %d, exiting early\n",
             data_vec_pos,end);
      exit(1);
   }

   double prior = 0.0;
   printf("computing mass prior...\n");
   for (int i = 0; i < start; ++i)
      prior += data_vec[i].thrust;
   prior /= start;

   printf("prior mass collected: %lf\n",prior);

   double post = 0.0;
   printf("computing mass post...\n");
   for (int i = end; i < data_vec_pos; ++i)
      post += data_vec[i].thrust;
   post /= (data_vec_pos - end);

   printf("post mass collected: %lf\n",post);

   double mB = (prior - post) / GRAVITY;
   printf("mass burned is: %lf\n",mB);

   struct data_sample * samples = &data_vec[start];
   int count = (end - start);
   double * M =  malloc(sizeof(double)*(count));
   double * best_M = malloc(sizeof(double)*(count));
   double best_Ve = 0.0;
   double best_err = 99999999.0;

   for (double Ve = 0.0; Ve < 3000.0; Ve += 0.00625) {

      compute_matching_M(samples,M,Ve,count);
      double err = E(samples,M,Ve,mB,count);

      if (err < best_err) {
         SWAP(M,best_M);
         best_Ve = Ve;
         best_err = err;
      }

   }

   printf("best_err: %e\n",best_err);

   FILE * prior_data = fopen("prior_data.txt","w");
   FILE * post_data = fopen("post_data.txt","w");

   fprintf(prior_data,"seconds\tthrust\n");
   fprintf(post_data,"seconds\tthrust\n");

   double time = 0.0;
   for (int i = 0; i < data_vec_pos; ++i) {

      fprintf(prior_data,"%.20lf\t%.20lf\n",
              time,data_vec[i].thrust);

      if (i < start || i >= end)
         fprintf(post_data,"%.20lf\t%.20lf\n",
                 time,data_vec[i].thrust);
      else
         fprintf(post_data,"%.20lf\t%.20lf\n",
                 time,best_M[i-start]*best_Ve);

      time += data_vec[i].sec;
   }

   fclose(prior_data);
   fclose(post_data);

   free(M);
   free(best_M);
   free(data_vec);
   return 0;
}
