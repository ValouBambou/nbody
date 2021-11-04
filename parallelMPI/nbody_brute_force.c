/*
** nbody_brute_force.c - nbody simulation using the brute-force algorithm (O(n*n))
**
**/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>
#include "nbody_tools.h"

#ifdef DISPLAY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#include "ui.h"
#include "nbody.h"
#include <mpi.h>

FILE *f_out = NULL;

int nparticles = 10; /* number of particles */
float T_FINAL = 1.0; /* simulation end time */
particle_t *particles;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

void init()
{
  /* Nothing to do */
}

#ifdef DISPLAY
Display *theDisplay; /* These three variables are required to open the */
GC theGC;            /* particle plotting window.  They are externally */
Window theMain;      /* declared in ui.h but are also required here.   */
#endif

/* compute the force that a particle with position (x_pos, y_pos) and mass 'mass'
 * applies to particle p
 */
void compute_force(particle_t *p, double x_pos, double y_pos, double mass)
{
  double x_sep, y_sep, dist_sq, grav_base;

  x_sep = x_pos - p->x_pos;
  y_sep = y_pos - p->y_pos;
  dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);

  /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
  grav_base = GRAV_CONSTANT * (p->mass) * (mass) / dist_sq;

  p->x_force += grav_base * x_sep;
  p->y_force += grav_base * y_sep;
}

/* compute the new position/velocity */
void move_particle(particle_t *p, double step)
{

  p->x_pos += (p->x_vel) * step;
  p->y_pos += (p->y_vel) * step;
  double x_acc = p->x_force / p->mass;
  double y_acc = p->y_force / p->mass;
  p->x_vel += x_acc * step;
  p->y_vel += y_acc * step;

  /* compute statistics */
  double cur_acc = (x_acc * x_acc + y_acc * y_acc);
  cur_acc = sqrt(cur_acc);
  double speed_sq = (p->x_vel) * (p->x_vel) + (p->y_vel) * (p->y_vel);
  double cur_speed = sqrt(speed_sq);

  sum_speed_sq += speed_sq;
  max_acc = MAX(max_acc, cur_acc);
  max_speed = MAX(max_speed, cur_speed);
}

double *xforces, *yforces;
int* local_sizes, *displs;
/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step, int lower_bound, int upper_bound, int np_local)
{
  /* First calculate force for particles. */
  double *xforces_local = malloc(sizeof(double) * np_local);
  double *yforces_local = malloc(sizeof(double) * np_local);
  int i;
  for (i = lower_bound; i < upper_bound; i++)
  {
    int j;
    particles[i].x_force = 0;
    particles[i].y_force = 0;
    for (j = 0; j < nparticles; j++)
    {
      particle_t *p = &particles[j];
      /* compute the force of particle j on particle i */
      compute_force(&particles[i], p->x_pos, p->y_pos, p->mass);
    }
    xforces_local[i - lower_bound] = particles[i].x_force;
    yforces_local[i - lower_bound] = particles[i].y_force;
  }
  // ask for forces computed by other processes
  MPI_Allgatherv(xforces_local, np_local, MPI_DOUBLE, xforces, local_sizes, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(yforces_local, np_local, MPI_DOUBLE, yforces, local_sizes, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  free(xforces_local);
  free(yforces_local);

  /* then move all particles and return statistics */
  for (i = 0; i < nparticles; i++)
  {
    particles[i].x_force = xforces[i];
    particles[i].y_force = yforces[i];
    move_particle(&particles[i], step);
  }
}

/* display all the particles */
void draw_all_particles()
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    int x = POS_TO_SCREEN(particles[i].x_pos);
    int y = POS_TO_SCREEN(particles[i].y_pos);
    draw_point(x, y);
  }
}

void print_all_particles(FILE *f)
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    particle_t *p = &particles[i];
    fprintf(f, "particle={pos=(%f,%f), vel=(%f,%f)}\n", p->x_pos, p->y_pos, p->x_vel, p->y_vel);
  }
}

void run_simulation(int rank, int size)
{
  int lower_bound, upper_bound;
  lower_bound = rank * (nparticles / size);
  upper_bound = (rank == size - 1) ? nparticles : lower_bound + (nparticles / size);
  int np_local = upper_bound - lower_bound; // nparticle for each MPI process
  double t = 0.0, dt = 0.01;
  MPI_Allgather(&np_local, 1, MPI_INT, local_sizes, 1, MPI_INT, MPI_COMM_WORLD); //   local_sizes[rank] = np_local; for each process
  MPI_Allgather(&lower_bound, 1, MPI_INT, displs, 1, MPI_INT, MPI_COMM_WORLD); //   displs[rank] = lower_bound; for each process

  while (t < T_FINAL && nparticles > 0)
  {
    /* Update time. */
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt, lower_bound, upper_bound, np_local);

    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */

    dt = 0.1 * max_speed / max_acc;

    /* Plot the movement of the particle */
#if DISPLAY
    clear_display();
    draw_all_particles();
    flush_display();
#endif
  }
}

/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char **argv)
{
  if (argc >= 2)
  {
    nparticles = atoi(argv[1]);
  }
  if (argc == 3)
  {
    T_FINAL = atof(argv[2]);
  }
  MPI_Init(&argc, &argv);

  init();

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Allocate global shared arrays for the particles data set. */
  particles = malloc(sizeof(particle_t) * nparticles);
  local_sizes = malloc(sizeof(int) * size);
  displs = malloc(sizeof(int) * size);
  xforces = malloc(sizeof(double) * nparticles);
  yforces = malloc(sizeof(double) * nparticles);
  all_init_particles(nparticles, particles);

  /* Initialize thread data structures */
#ifdef DISPLAY
  /* Open an X window to display the particles */
  simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

  struct timeval t1, t2;
  gettimeofday(&t1, NULL);

  /* Main thread starts simulation ... */

  run_simulation(rank, size);

  gettimeofday(&t2, NULL);

  double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

#ifdef DUMP_RESULT
  FILE *f_out = fopen("particles.log", "w");
  assert(f_out);
  print_all_particles(f_out);
  fclose(f_out);
#endif
  
  if (rank == 0)
  {
    printf("-----------------------------\n");
    printf("nparticles: %d\n", nparticles);
    printf("T_FINAL: %f\n", T_FINAL);
    printf("-----------------------------\n");
    printf("Simulation took %lf s to complete\n", duration);
  }

  free(xforces);
  free(yforces);
  free(local_sizes);

#ifdef DISPLAY
  clear_display();
  draw_all_particles();
  flush_display();

  printf("Hit return to close the window.");

  getchar();
  /* Close the X window used to display the particles */
  XCloseDisplay(theDisplay);

#endif
  free(particles);
  MPI_Finalize();

  return 0;
}
