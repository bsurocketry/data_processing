
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

typedef struct constant_values {
   double A_star;
   double Pe;
   double Po;
   double Pa;
   double Ae;
   double k;
   double R;
   double pp;
} constant_values;

static void compute_matching_M(struct data_sample * samples, double * M, double Ve, int count) {
   double accum = 0.0; /* previous mass burned */

   for (int i = 0; i < count; ++i) {

      M[i] = (samples[i].thrust + (accum*GRAVITY)) / Ve;

      accum += M[i] * samples[i].sec;

   }
}

static void compute_matching_Po(double * M, double * Po, double To, constant_values * c, int count) {

   for (int i = 0; i < count; ++i) {
      Po[i] = M[i] / (c->A_star*sqrt((c->k/(c->R*To))*(pow(2.0/(c->k+1),(c->k+1.0)/(2*(c->k-1))))));
   }

}

static double Po_E(struct data_sample * samples, double * Po, constant_values * c, int count) {

   double err = 0.0;

   for (int i = 0; i < count; ++i) {

      double lhs = pow((samples[i].thrust-(c->Pe-c->Pa)*c->Ae)/c->A_star,2.0);

      double rhs = ((2.0*c->k*c->k)/(c->k-1))*pow(2.0/(c->k+1),(c->k+1)/(c->k-1));

      rhs *= (Po[i]*Po[i])-(pow(c->Pe,(c->k-1)/c->k)*pow(c->Po,2.0-(c->k-1)/c->k));

      /* squared error */
      err += (rhs - lhs) * (rhs - lhs);

   }

   return err;
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

#define MIN(a,b) \
   ((a) < (b) ? (a) : (b))

static inline double triangle(double A, double fA, double B, double fB) {
   double height = fB - fA;
   double width  = (B - A);
   double rect = width * MIN(fA,fB);
   return ((height * width) / 2.0) + rect;
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

   double prior_time = 0.0;
   for (int i = 0; i < start; ++i)
      prior_time += data_vec[i].sec;

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
   printf("best_Ve:  %lf\n",best_Ve);

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
      else {
         fprintf(post_data,"%.20lf\t%.20lf\n",
                 time,best_M[i-start]*best_Ve);

         /* rewrite history! xD */
         data_vec[i].thrust = best_M[i-start]*best_Ve;
      }

      time += data_vec[i].sec;
   }

   double total_impulse = 0.0;
   time = 0.0;
   for (int i = 0; i < count-1; ++i) {
      time += samples[i].sec;
      total_impulse += triangle(time,
                                best_M[i]*best_Ve,
                                time + samples[i+1].sec,
                                best_M[i+1]*best_Ve);
   }

   printf("total impulse: %lf\n",total_impulse);

   /* now we will begin the process of solving for Po and k
    * - these are pressure of exhaust and burn rate.        */

   constant_values c;
   bzero(&c,sizeof(constant_values));
   c.Pe = 1.013 * 1e5;
   c.k = 1.21;
   c.Pa = c.Pe;
   c.Ae = 0.0;
   c.A_star = (M_PI / 4.0) * 10e-5;
   c.R = 208.37;
   c.pp = 1820.0;

   double best_To = 0.0;
   best_err = 99999999999999999999999999999999.0;
   double * Po = malloc(sizeof(double)*(count));
   double * best_Po = malloc(sizeof(double)*(count));
   for (double To = 0.0; To < 10000.0; To += 0.00625) {
      compute_matching_Po(best_M,Po,To,&c,count);
      double err = Po_E(samples,Po,&c,count);

      if (err < best_err) {
         SWAP(Po,best_Po);
         best_To = To;
         best_err = err;
      }
   }

   printf("best_err: %e\n",best_err);
   printf("best_To: %lf\n",best_To);

   FILE * Po_data = fopen("chamber_pressure.txt","w");
   fprintf(Po_data,"seconds\tpressure(Pa)\n");

   time = prior_time;
   for (int i = 0; i < count; ++i) {

      fprintf(Po_data,"%.20lf\t%.20lf\n",
              time,best_Po[i]);

      time += samples[i].sec;
   }

   fclose(prior_data);
   fclose(post_data);
   fclose(Po_data);

   free(M);
   free(best_M);
   free(data_vec);
   free(Po);
   free(best_Po);
   return 0;
}
