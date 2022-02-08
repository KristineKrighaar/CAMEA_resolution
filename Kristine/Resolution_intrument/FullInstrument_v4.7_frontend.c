/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: FullInstrument_v4.7_frontend.instr (PSI_CAMEA)
 * Date:       Wed Oct 27 13:42:07 2021
 * File:       ./FullInstrument_v4.7_frontend.c
 * Compile:    cc -o PSI_CAMEA.out ./FullInstrument_v4.7_frontend.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.7.1 - Oct. 04, 2021"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 2.7.1
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */

#ifndef WIN32
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#endif

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.7.1 - Oct. 04, 2021"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Oct. 04, 2021"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.7.1"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */

/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
/* When using c99 in the CFLAGS, some of these consts
   are lost... Perhaps we should in fact include everything from
   https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
*/
#  define PI 3.14159265358979323846
#  define M_PI PI
#  define M_PI_2 M_PI/2.0
#  define M_PI_4 M_PI/4.0
#  define M_1_PI 1.0/M_PI
#  define M_2_PI 2*M_1_PI
#  define M_2_SQRTPI 2/sqrt(M_PI)
#  define M_SQRT2 sqrt(2)
#  define M_SQRT1_2 sqrt(1/2)
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("\nINSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 712 "./FullInstrument_v4.7_frontend.c"

#line 1 "mcstas-r.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundant with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; goto mcabsorbComp; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlz=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; } \
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlx=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; } \
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnly=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; } \
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 945 "./FullInstrument_v4.7_frontend.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#include <sys/stat.h>

#ifdef _WIN32 
#include <direct.h>
# define  mkdir( D, M )   _mkdir( D ) 
#endif 

#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Reduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;
  
  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);
      
    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+"); 
    
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n", 
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m)
    return(detector);
  
  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];
  
  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m) 
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin); 
      else x = 0;
      if (detector.n && detector.p) 
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin); 
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw) 
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */
      
      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }
      
      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
    
    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH, 
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }
  
  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  
  return(detector);
  
} /* mcdetector_statistics */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre, 
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out_backend: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii, mcinfo(stdout)
*******************************************************************************/
static void mcruninfo_out_backend(char *pre, FILE *f, int info)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
      if (!info){
          (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
          fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
      }else{
        /*if an info run, some variables might not have values. Flag these by "NULL"*/
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
            /* ... those with defautl values*/
            (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
            fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        }else{
            /* ... and those without */
            fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
} /* mcruninfo_out_backend */

/************************
* wrapper function to mcruninfo_out_backend
*  Regular runs use this whereas the single call from mcinfo is directly to the backend
*************************/
static void mcruninfo_out(char *pre, FILE *f){
    mcruninfo_out_backend(pre,f,0);
}

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;
  
  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" : 
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f, 
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" : 
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre, 
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
    
  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  /* Write data set information to simulation description file. */
  MPI_MASTER(
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n", 
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);
      
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
  
}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
      }
      fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2, 
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0, 
          outfile, detector.istransposed);
      }
      fclose(outfile);
      
      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",  
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=CHAR_BUF_LENGTH; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;
  
  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);
  
  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);
  
    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }
  
  return(ret);
  
} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: mcsiminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];
  
  if (!f || mcdisable_output_files) return;
  
  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */ 
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */
      
      NXclosegroup(f); /* simulation */
    } /* NXsimulation */
    
    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];
  if (!f || !detector.m || mcdisable_output_files) return;
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
    
    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
    
      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" : 
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" : 
                 "xylimits", detector.limits);
      nxprintattr(f, "variables", 
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");
        
      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */
  
} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[CHAR_BUF_LENGTH];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);
    
    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{
  
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;
  
  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);
  
  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;
  
  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */
  
  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",  
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }
  
  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];
  
  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
  
      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar, 
          1, detector.m, detector.xmin, detector.xmax);
          
        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar, 
          2, detector.n, detector.ymin, detector.ymax);
          
        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar, 
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */
      
      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);
      
      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may. 
*   Then the number of recv/send is not constant along nodes, and simulation stalls.  
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );
  
  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
	    
	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );   
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );
  
  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  
#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
{
  int exists=0;
  int index;
  
  /* check format */      
  if (!mcformat || !strlen(mcformat) 
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE") 
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }
  
  /* open the SIM file if not defined yet */
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));
    
} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));
  
} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];
  
  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[strcspn(xvar,"\n\r ")]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[strcspn(yvar,"\n\r ")]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));
  
} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;
  
  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,labs(m),1,labs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);
  
  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out_backend("  ", stdout,1);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*M_PI/24),y,z+r*cos(i*2*M_PI/24),
                    x+r*sin((i+1)*2*M_PI/24),y,z+r*cos((i+1)*2*M_PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*M_PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, M_PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_invert: Matrix inversion, in case a Rotatoin is used to represent a
* a general non-orthonormal matrix.
*******************************************************************************/
void
rot_invert(Rotation t1, Rotation t2)
{
  Rotation cofactors;
  int r,c;
  double det=0;
  for (r=0;r<3;r++){
    for (c=0;c<3;c++){
      /*this algorithm automatically takes care of the sign changes in computing cofactors*/
      cofactors[r][c]=t1[(r+1) % 3][(c+1) % 3]*t1[(r+2) % 3][(c+2) % 3] - t1[(r+2) % 3][(c+1) % 3]*t1[(r+1) % 3][(c +2) % 3] ;
    }
  }
  det=t1[0][0]*cofactors[0][0] + t1[0][1]*cofactors[0][1] + t1[0][2]*cofactors[0][2];
  if(det==0){
    fprintf(stderr,"Warning: matrix not invertable\n");
  }
  rot_transpose(cofactors,t2);

  /*the adjoint matrix should now be scaled by 1/det to get the inverse*/
  for (r=0;r<3;r++){
    for (c=0;c<3;c++){
      t2[r][c]=t2[r][c]/det;
    }
  }
}
/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order_old(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

int solve_2nd_order(double *t0, double *t1, double A, double B, double C){
  int retval=0;
  double sign=copysign(1.0,B);
  double dt0,dt1;

  dt0=0;
  dt1=0;
  *t0;
  if(t1){ *t1=0;}

  /*protect against rounding errors by locally equating DBL_EPSILON with 0*/
  if (fabs(A)<DBL_EPSILON){
    A=0;
  }
  if (fabs(B)<DBL_EPSILON){
    B=0;
  }
  if (fabs(C)<DBL_EPSILON){
    C=0;
  }

  /*check if coefficient are sane*/
  if( A==0  && B==0){
    retval=0;
  }else{
    if(A==0){
      /*equation is linear*/
      dt0=-C/B;
      retval=1;
    }else if (C==0){
      /*one root is 0*/
      if(sign<0){
        dt0=0;dt1=-B/A;
      }else{
        dt0=-B/A;dt1=0;
      }
      retval=2;
    }else{
      /*a regular 2nd order eq. Also works out fine for B==0.*/
      double D;
      D=B*B-4*A*C;
      if (D>=0){
        dt0=(-B - sign*sqrt(B*B-4*A*C))/(2*A);
        dt1=C/(A*dt0);
        retval=2;
      }else{
        /*no real roots*/
        retval=0;
      }
    }
    /*sort the solutions*/
    if (retval==1){
      /*put both solutions in t0 and t1*/
      *t0=dt0;
      if(t1) *t1=dt1;
    }else{
      /*we have two solutions*/
      /*swap if both are positive and t1 smaller than t0 or t1 the only positive*/
      int swap=0;
      if(dt1>0 && ( dt1<dt0 || dt0<=0) ){
        swap=1;
      }
      if (swap){
        *t0=dt1;
        if(t1) *t1=dt0;
      }else{
        *t0=dt0;
        if(t1) *t1=dt0;
      }
    }

  }
  return retval;

} /*solve_2nd_order_improved*/


/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  
#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif
  srandom(mcseed);

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 5083 "./FullInstrument_v4.7_frontend.c"

#line 1 "mcstas-r.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999. 
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98 
  *******************************************************************************/
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1 
 *******************************************************************************/
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int
plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */

#line 5443 "./FullInstrument_v4.7_frontend.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "PSI_CAMEA";
char mcinstrument_source[] = "FullInstrument_v4.7_frontend.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Source_gen4'. */
#line 132 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

/*maximum number of rows to rebin a table = 1M*/
enum { mcread_table_rebin_maxsize = 1000000 };

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 5052 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value should just be used as
     * if the table had been read from disk. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we should proceed with freeing the memory
     * associated with the table item - otherwise only decrement the reference counter since there are more references
     * that may need it.*/

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref = ((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found and no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found and the reference counter is 1.
                         * This means we should garbage collect. Move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            /* item not found, and so should be garbage collected. This could be the case if freeing a
             * Table that has been constructed from code - not read from file. Return 0x1 to flag it for
             * collection.*/
            return (void *) 0x1 ;
    }
}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    return Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && mcinstrument_source[0] != '\0' && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe[0] != '\0' && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        // printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array*1.5;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename[0] != '\0' ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      /*return early if the rebinned table will become too large*/
      if (Length_Table > mcread_table_rebin_maxsize){
        fprintf(stderr,"WARNING: (Table_Rebin): Rebinning table from %s would exceed 1M rows. Skipping.\n", Table->filename); 
        return(Table->rows*Table->columns);
      }
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index-1,0);
    X2 = Table_Index(Table,Index  ,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index <= Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1, j);
  Y2 = Table_Index(Table,Index  , j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table. First Call Table_File_list_gc. If this returns
*   non-zero it means there are more refernces to the table, and so the table
*   should not bee freed.
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename[0] != '\0' ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* first allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /*if the block is empty - don't store it*/
      if (nelements>0){
          /* if t_Table array is not long enough, expand and realocate */
          if (block_number >= allocated-1) {
              allocated += 256;
              Table_Array = (t_Table *)realloc(Table_Array,
                      allocated*sizeof(t_Table));
              if (!Table_Array) {
                  fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
                          allocated*sizeof(t_Table));
                  *blocks = 0;
                  return (NULL);
              }
          }
          /* store it into t_Table array */
          //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
          Table_Array[block_number-1] = Table;
      }
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index;
    if (!Table) return;
    for (index=0;index < Table[0].array_length; index++){
            Table_Free(&Table[index]);
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


#ifndef SOURCE_GEN_DEF
#define SOURCE_GEN_DEF
/*******************************************************************************
* str_dup_numeric: makes a clean copy of a string and allocate as numeric
*******************************************************************************/
char *str_dup_numeric(char *orig)
{
  long i;
  char *valid;

  if (!orig || !strlen(orig)) return(NULL);

  for (i=0; i < strlen(orig); i++)
  {
    if ( (orig[i] > 122)
      || (orig[i] < 32)
      || (strchr("!\"#$%&'()*,:;<=>?@[\\]^`/ ", orig[i]) != NULL) )
    {
      orig[i] = ' ';
    }
  }
  orig[i] = '\0';
  /* now skip spaces */
  for (i=0; i < strlen(orig); i++) {
    if (*orig == ' ') orig++;
    else break;
  }

  return(orig);
} /* str_dup_numeric */
#endif
#line 6938 "./FullInstrument_v4.7_frontend.c"

/* Shared user declarations for all components 'Elliptic_guide_gravity'. */
#line 187 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Depends on read_table-lib
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/


#ifndef REF_LIB_H
#define REF_LIB_H "$Revision$"

void StdReflecFunc(double, double*, double*);
void TableReflecFunc(double, t_Table*, double*);

#endif

/* end of ref-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.c
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Variable names have prefix 'mc_ref_' for 'McStas Reflection' 
* to avoid conflicts
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/

#ifndef REF_LIB_H
#include "ref-lib.h"
#endif

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#include "read_table-lib.c"
#endif

/****************************************************************************
* void StdReflecFunc(double q, double *par, double *r)
* 
* The McStas standard analytic parametrization of the reflectivity.
* The parameters are:
* R0:      [1]    Low-angle reflectivity
* Qc:      [AA-1] Critical scattering vector
* alpha:   [AA]   Slope of reflectivity
* m:       [1]    m-value of material. Zero means completely absorbing.
* W:       [AA-1] Width of supermirror cut-off
*****************************************************************************/
void StdReflecFunc(double mc_pol_q, double *mc_pol_par, double *mc_pol_r) {
    double R0    = mc_pol_par[0];
    double Qc    = mc_pol_par[1];
    double alpha = mc_pol_par[2];
    double m     = mc_pol_par[3];
    double W     = mc_pol_par[4];
    double beta  = 0;
    mc_pol_q     = fabs(mc_pol_q);
    double arg;
        
    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       double m_value=m*0.9853+0.1978;
       double W=-0.0002*m_value+0.0022;
       double alpha=0.2304*m_value+5.0944;
       double beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));
    */  
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	alpha=m;
	beta=0;
      }
    }
    
    arg = W > 0 ? (mc_pol_q - m*Qc)/W : 11;

    if (arg > 10 || m <= 0 || Qc <=0 || R0 <= 0) {
      *mc_pol_r = 0;
      return;
    }
    
    if (m < 1) { Qc *= m; m=1; }
    
    if(mc_pol_q <= Qc) {      
      *mc_pol_r = R0;
      return;
    }
    
    
    *mc_pol_r = R0*0.5*(1 - tanh(arg))*(1 - alpha*(mc_pol_q - Qc) + beta*(mc_pol_q - Qc)*(mc_pol_q - Qc));
    
    return;
  }

/****************************************************************************
* void TableReflecFunc(double q, t_Table *par, double *r) {
* 
* Looks up the reflectivity in a table using the routines in read_table-lib.
*****************************************************************************/
void TableReflecFunc(double mc_pol_q, t_Table *mc_pol_par, double *mc_pol_r) {
    
  *mc_pol_r = Table_Value(*mc_pol_par, mc_pol_q, 1);
  if(*mc_pol_r>1)
    *mc_pol_r = 1;
  return;
}

/* end of ref-lib.c */


///////////////////////////////////////////////////////////////////////////
/////////////// local structs and enums
///////////////////////////////////////////////////////////////////////////

/**
Sides of the guide
*/
enum Side {RightSide,TopSide,LeftSide,BottomSide,None};

/**
The type of the collision is set in the collision function
and decide the functions called in trace()
	Reflex (TODO change this name) calls the reflection function
	Absorb calls the built in ABSORB funtion.
	LeaveGuide calls break and end the calculations in this component
	EnterGuide does nothing
*/
enum CollisionType {Reflex,Absorb,LeaveGuide,EnterGuide};

/**
	The Mirror type sets the CollisionType of particles colliding on the mirror
*/
enum MirrorType {MirrorTypeReflection,MirrorTypeTransparent,MirrorTypeabsorption};

// enum IntersectionType {Reflex,Absorb,Transparent,Leave,Enter};

/**
	Collision between guide and the particle
	contain infomation on the time to the next collision,
	which side of the guide it is on and whether this part of the guide
	is a perfect or approximated ellipse.
*/
struct Intersection
{
	double delta_time_to_next_collision; 
	enum Side side;	// A number from 0 to 4 (4 being an error warning)
	int ApproxOn;
	enum CollisionType collisionType;
};

/**
	Static Guide information (SGI)
	contain information on the guide, the ellipses and the mirrors on all sides
*/
struct SGI
{
	// guide infomation
	double Length;
	double entranceHorizontalWidth, entranceVerticalWidth;	
	double exitHorizontalWidth, exitVerticalWidth;	

	// ellipses infomation
	double ellipseMajorAxis[4],		ellipseMinorAxis[4];
	double ellipseMajorOffset[4],	ellipseMinorOffset[4];

	// mirror infomation
	double R0Arr[4]; 
	double QcArr[4]; 
	double alphaArr[4]; 
	double mArr[4]; 
	double WArr[4];

	// mirror type
	enum MirrorType InnerSide[4];
	enum MirrorType OuterSide[4];

	// selene
	int EnclosingBoxOn;
    double xArray[8];
    double yArray[8];
    double zArray[8];

    // segmentation
    int numberOfSegments;
    int enableSegments;
    double *mValuesright;
	double *mValuesleft;
	double *mValuestop;
	double *mValuesbottom;
    double *segLength;

    int verboseSetting;
};


///////////////////////////////////////////////////////////////////////////
/////////////// Error Handling Functions
///////////////////////////////////////////////////////////////////////////

/**
	If a user input is less than zero and hence doesn't allow for a well
	define geomtric of the guide or physical values for mirrors
	@param var is the input varible there the error occurred [text]
*/
int guide_elliptical_illegalInputLessThanZero(char var[],int verbose){
	if (verbose)
		printf("The user defined variable %s in %s has an illegal value"
				" less than zero\n",var,"Elliptic_guide_gravity");
	return 1;
}

/**
	The first focal point is in and the second is out.
	If -in-out > L then they would change position as the 
	first and second focal points. This is 
	@param in,out is the input varible there the error occurred [text]
*/
int guide_elliptical_illegalInputFocalPointsHyperbola(
			char in[],char out[],
			double inValue,double outValue, int verbose){
	if (verbose){
		printf("The user defined length of the guide, length \
				and the focal points %s and %s does not result \
				in an well defined ellipse. swap the focal points \
				or increase L, %s or %s to fix this problem\n",
				in,out,in,out);
		printf("The mininum length of the should be around %e\n",
				inValue+outValue+0.000001);
	}
	return 1;
}

/**
	Gives a warning if a part of the code is called that 
	should not be accessible if the algoritmes are working correctly
	Most likely errors are floating points and ill-defined cases
*/
void guide_elliptical_callCriticalWarning(char func[],int verbose){
	if (verbose)
		printf("A CRITICAL WARNING has been called inside %s by function %s." 
			"This is most likely due to a programming error \
			inside the component. \n",
			"Elliptic_guide_gravity",func);
	}

///////////////////////////////////////////////////////////////////////////
/////////////// Collision handling functions
///////////////////////////////////////////////////////////////////////////

int guide_elliptical_getMirrorTypeFromInput(char * input,int verbose){
	int type = -1;
	char* r1 = "reflection"; char* r2 = "reflect"; char* r3 = "r";
	char* a1 = "absorption"; char* a2 = "absorb"; char* a3 = "a";
	char* t1 = "transparant";char* t2 = "trans"; char* t3 = "t";
	if (strcmp (input, r1) == 0
	||  strcmp (input, r2) == 0 
	||  strcmp (input, r3) == 0)
		type = MirrorTypeReflection;
	if (strcmp (input, a1) == 0
	||  strcmp (input, a2) == 0 
	||  strcmp (input, a3) == 0)
		type = MirrorTypeabsorption;
	if (strcmp (input, t1) == 0
	||  strcmp (input, t2) == 0 
	||  strcmp (input, t3) == 0)
		type = MirrorTypeTransparent;
	if ( type == -1 && verbose)
		printf( "Following string is not a valid type of a mirror: %s,"
				"use reflection,absorption or transparant. \n" ,input);

	return type;
	}

///////////////////////////////////////////////////////////////////////////
/////////////// Collision functions
///////////////////////////////////////////////////////////////////////////	

/**
	Find the intersection between the neutron and the ellipse using newton method.
	As there is up to 4 solution to this problem, and only the 
	smallest positive root is the physical solution. Using the tuning points 
	it is possible to look the only the potential roots to speed up calculations.

	@param coef; A pointer to the array holding the coeffecients 
			for the 4th order polynomial.
	@param startPosition, The default starting point for newton method. [s]
	@param limit; A point after all the roots of the polynial. [s]
	@param solution A pointer which will hold the physical solution 
			if this function return true.
	@return; return 1 if the physical solution is found. [boolean] 
*/

double guide_elliptical_foverdf(double *coefficients,double currentPoint){
	double numerator= coefficients[0]*currentPoint*currentPoint*currentPoint*currentPoint
					+ coefficients[1]*currentPoint*currentPoint*currentPoint
					+ coefficients[2]*currentPoint*currentPoint
					+ coefficients[3]*currentPoint
					+ coefficients[4];
	double denominator=4*coefficients[0]*currentPoint*currentPoint*currentPoint
					+ 3*coefficients[1]*currentPoint*currentPoint
					+ 2*coefficients[2]*currentPoint
					+ coefficients[3];
	return numerator/denominator;
}

int guide_elliptical_newtonRapsonsMethod4thOrder(
		double *coefficients,double *solution,double startingPoint,
		double tolerance,double max_iterations){

	double numerator;
	double denominator;
	double t_previous;
	double t = startingPoint;
	int iteration = 0;

	do {
		t_previous = t;
		t = t_previous - guide_elliptical_foverdf(coefficients,t);
		iteration++;
	} while( fabs(t-t_previous) > tolerance && iteration < max_iterations );
	if( iteration == max_iterations ) { return 0; }
	else 			{ *solution = t;	return 1; }
}



int guide_elliptical_findNeutronEllipseIntersection(
				double *coef,double startPosition,
				double limit,double *solution){

	// in the case of no gravity
	if(coef[0] == 0 & coef[1] == 0){
		double t1=0;
		double t2=0;
		int boolean = solve_2nd_order(&t1,&t2,coef[2],coef[3],coef[4]);

		if ( t1 > startPosition ){ *solution = t1; }
		if ( t2 > startPosition ){ *solution = t2; }
		return boolean;
	}

	double tol = 1e-15;
	double max_iter = 1e3;
	double turningP1,turningP2;
	
	double sp = startPosition;
	int inside;
	if ( coef[0]*sp*sp*sp*sp
		+coef[1]*sp*sp*sp
		+coef[2]*sp*sp
		+coef[3]*sp
		+coef[4] < 0)
		 inside = 1;
	else inside = 0;

	int boolean = solve_2nd_order(
						&turningP1,&turningP2,
						12*coef[0],6*coef[1],2*coef[2]);

	double t1=0,t2=0;
	double ss=100;

	if( inside ){
		if(boolean) guide_elliptical_newtonRapsonsMethod4thOrder(coef,&t1,turningP1,tol,max_iter);
		guide_elliptical_newtonRapsonsMethod4thOrder(coef,&t2,limit,tol,max_iter);
	}
	else{
		if(boolean) guide_elliptical_newtonRapsonsMethod4thOrder(coef,&t1,turningP2,tol,max_iter);
		guide_elliptical_newtonRapsonsMethod4thOrder(coef,&t2,startPosition,tol,max_iter);
	}

	if (ss > t1 && t1 > 1e-15)	ss = t1;
	if (ss > t2 && t2 > 1e-15)	ss = t2;
	*solution = ss;
	
	return 1;
}


int guide_elliptical_handleGuideIntersection(
		double x, double y, double z,
		double vx,double vy,double vz,
		double Gx,double Gy,double Gz,
		struct SGI *guideInfo,
		struct Intersection *currentCollision){
	//
	double horExS = 1/( guideInfo->ellipseMinorAxis[RightSide]
					   *guideInfo->ellipseMinorAxis[RightSide]);
	double horEzS = 1/( guideInfo->ellipseMajorAxis[RightSide]
				 	   *guideInfo->ellipseMajorAxis[RightSide]);
	double hordiffx = x-guideInfo->ellipseMinorOffset[RightSide]; 
	double hordiffz = z-guideInfo->ellipseMajorOffset[RightSide]; 

	double horAlpha = ( Gx*Gx*horExS + Gz*Gz*horEzS )/4;
	double horBeta  = ( Gx*vx*horExS + Gz*vz*horEzS );
	double horGamma = horExS*vx*vx + horEzS*vz*vz 
					 + horExS*Gx*hordiffx + horEzS*Gz*hordiffz; 
	double horDelta = 2*horExS*vx*hordiffx + 2*horEzS*vz*hordiffz;
	double horEpsilon = horExS*hordiffx*hordiffx + horEzS*hordiffz*hordiffz - 1; 

	double horCoefficients[5] = {horAlpha,horBeta,horGamma,horDelta,horEpsilon};

	double verEyS = 1/( guideInfo->ellipseMinorAxis[TopSide]
					   *guideInfo->ellipseMinorAxis[TopSide]);
	double verEzS = 1/( guideInfo->ellipseMajorAxis[TopSide]
				 	   *guideInfo->ellipseMajorAxis[TopSide]);
	double verdiffy = y-guideInfo->ellipseMinorOffset[TopSide];
	double verdiffz = z-guideInfo->ellipseMajorOffset[TopSide];

	double verAlpha = ( Gy*Gy*verEyS + Gz*Gz*verEzS )/4;
	double verBeta  = ( Gy*vy*verEyS + Gz*vz*verEzS );
	double verGamma = verEyS*vy*vy + verEzS*vz*vz 
					+ verEyS*Gy*verdiffy + verEzS*Gz*verdiffz;
	double verDelta = 2*verEyS*vy*verdiffy + 2*verEzS*vz*verdiffz;
	double verEpsilon = verEyS*verdiffy*verdiffy + verEzS*verdiffz*verdiffz - 1;

	double verCoefficients[5] = {verAlpha,verBeta,verGamma,verDelta,verEpsilon};


	double upperlimit;
	double startingPoint = 1e-15;

	int boolean;
	// Horizontal
	double solutionH = 0;
	solve_2nd_order(
			&upperlimit,NULL,
			-0.5*Gz,-vz,2*guideInfo->ellipseMajorAxis[RightSide]-z);
	int booleanH = guide_elliptical_findNeutronEllipseIntersection(
						horCoefficients,startingPoint,upperlimit,&solutionH);
	// Vertical
	double solutionV = 0;
	solve_2nd_order(
			&upperlimit,NULL,
			-0.5*Gz,-vz,2*guideInfo->ellipseMajorAxis[TopSide]-z);
	int booleanV = guide_elliptical_findNeutronEllipseIntersection(
						verCoefficients,startingPoint,upperlimit,&solutionV);

	if (solutionH <= 0)
			currentCollision->delta_time_to_next_collision = solutionV;
	else if (solutionV <= 0)
			currentCollision->delta_time_to_next_collision = solutionH;
	else if (fabs(solutionH - solutionV) < 1e-12) return 0;
	else if (solutionH < solutionV){
		currentCollision->delta_time_to_next_collision = solutionH; 
		boolean = booleanH; 
		}
	else{
		currentCollision->delta_time_to_next_collision = solutionV; 
		boolean = booleanV; 
		}

	double tside = currentCollision->delta_time_to_next_collision;
	double xside = x + vx*tside + 0.5*Gx*tside*tside;
	double yside = y + vy*tside + 0.5*Gy*tside*tside;
	double zside = z + vz*tside + 0.5*Gz*tside*tside;

	double xfactor =
		2*sqrt(1 - ( (zside-guideInfo->ellipseMajorOffset[RightSide])
					*(zside-guideInfo->ellipseMajorOffset[RightSide]) 
					)/(guideInfo->ellipseMajorAxis[RightSide] 
					*guideInfo->ellipseMajorAxis[RightSide] ) 
		)*guideInfo->ellipseMinorAxis[RightSide];

	double yfactor =
		2*sqrt(1 - ( (zside-guideInfo->ellipseMajorOffset[BottomSide])
					*(zside-guideInfo->ellipseMajorOffset[BottomSide]) 
					)/(guideInfo->ellipseMajorAxis[BottomSide] 
					*guideInfo->ellipseMajorAxis[BottomSide] )
		)*guideInfo->ellipseMinorAxis[BottomSide];

	xside = xside/xfactor;
	yside = yside/yfactor;
	if( fabs(yside) >= fabs(xside) ){
		if(y > 0)	currentCollision->side = TopSide;
		else 		currentCollision->side = BottomSide;
	}
	else{
		if(x < 0)	currentCollision->side = RightSide;
		else 		currentCollision->side = LeftSide;
	}
	if (tside < 1e-15) printf("low time is: %e\n",tside);

	return boolean;
}

/**
	Check if the neutron is within the guide using the sign 
	of the crossproduct between the two points,
	on each of the enclosing box surface and neutrons position.    
	
	@param x,y,z; position of the neutron. [m]
	@param guideInfo; pointer to the guide infomation holding structure.
	@return; return 1 if the neutron is inside the guide [boolean] 
*/

/*
int guide_elliptical_InsideEnclosingBox(double x,double y,double z,struct SGI *guideInfo){
	int guide_elliptical_IsPointInVolume(
			double *x,double *y,double *z,
			double px,double py,double pz){
		int guide_elliptical_WhichSide( 	double p1x,double p1y,double p1z,
						double p2x,double p2y,double p2z,
						double p3x,double p3y,double p3z,
						double px ,double py ,double pz ){
			
			double v1x = p1x - p2x, v1y = p1y-p2y, v1z = p1z-p2z;
			double v2x = p3x - p2x, v2y = p3y-p2y, v2z = p3z-p2z;
			double v3x = v2y*v1z-v2z*v1y; 
			double v3y = v2z*v1x-v2x*v1z;
			double v3z = v2x*v1y-v2y*v1x;

			return 0 >= v3x*(px-p1x)+v3y*(py-p1y)+v3z*(pz-p1z); 
		}

		if(	//front
			guide_elliptical_WhichSide(x[3],y[3],z[3],x[2],y[2],z[2],x[1],y[1],z[1],px,py,pz) &&
			guide_elliptical_WhichSide(x[1],y[1],z[1],x[0],y[0],z[0],x[3],y[3],z[3],px,py,pz) &&
			//back
			guide_elliptical_WhichSide(x[5],y[5],z[5],x[6],y[6],z[6],x[7],y[7],z[7],px,py,pz) &&
			guide_elliptical_WhichSide(x[7],y[7],z[7],x[4],y[4],z[4],x[5],y[5],z[5],px,py,pz) &&
			//right
			guide_elliptical_WhichSide(x[7],y[7],z[7],x[3],y[3],z[3],x[0],y[0],z[0],px,py,pz) &&
			guide_elliptical_WhichSide(x[0],y[0],z[0],x[4],y[4],z[4],x[7],y[7],z[7],px,py,pz) &&
			//left
			guide_elliptical_WhichSide(x[1],y[1],z[1],x[2],y[2],z[2],x[6],y[6],z[6],px,py,pz) &&
			guide_elliptical_WhichSide(x[6],y[6],z[6],x[5],y[5],z[5],x[1],y[1],z[1],px,py,pz) &&
			//top
			guide_elliptical_WhichSide(x[0],y[0],z[0],x[1],y[1],z[1],x[5],y[5],z[5],px,py,pz) &&
			guide_elliptical_WhichSide(x[5],y[5],z[5],x[4],y[4],z[4],x[0],y[0],z[0],px,py,pz) &&
			//bottom
			guide_elliptical_WhichSide(x[6],y[6],z[6],x[2],y[2],z[2],x[3],y[3],z[3],px,py,pz) &&
			guide_elliptical_WhichSide(x[3],y[3],z[3],x[7],y[7],z[7],x[6],y[6],z[6],px,py,pz) )
			 	return 1;
		else 	return 0;
	}
	return guide_elliptical_IsPointInVolume(
				guideInfo->xArray,guideInfo->yArray,guideInfo->zArray,x,y,z);
}
*/


///////////////////////////////////////////////////////////////////////////
/////////////// reflection functions
///////////////////////////////////////////////////////////////////////////
		

/**
	Calculate the new velocity vector for the particle colliding on 
	the inner side of the elliptic mirror and returns the loss-factor (TODO)
	
	@param pos_V0,pos_W0 Is the 2d position vector of the particle, 
			assumed to be a point on the ellipse. [m]
	@param pvel_V0,pvel_W0 Is the 2d velocity vector of the particle. [m/s]
	@param ellipse_V_axis_squared,ellipse_W_axis_squared
			are the axes of the ellipse. [m]
	@param ellipse_V_offset,ellipse_W_offset Is the 2d vector difference 
			between the ellipse coordinate system (center of the ellipse) 
			and the guide coordinate system [m]
	@param R0, Mvalue, Qc, W, Alpha #TODO 
		slaa beskrivelse af disse variabler i andre dokumenter 
		og hold dig til standarden.
	@return the new wieght of the package
*/
double guide_elliptical_ReflectionOnEllipticSurface(
			double pos_V,double pos_W,
			double *pvel_V,double *pvel_W,
			double ellipse_V_axis,double ellipse_W_axis,
			double ellipse_V_offset,double ellipse_W_offset,
			double R0, double Qc, double alpha, double Mvalue, double W)
{

	// Turns the velocity vector (vel_V0,vel_W0) into a local value
	double vel_V = *pvel_V;
	double vel_W = *pvel_W;

	// Galilean transformation of the particles start position 
	// to the ellipse coordinate system
	pos_V=pos_V-ellipse_V_offset;
	pos_W=pos_W-ellipse_W_offset;

	/*  
	*   If we reflect the velocity vector in the normal 
	*	to the ellipse in the point of intersection
	*   The resulting vector will be -f2, do to conservation of momentum.
	*   this result in the following equation
	*   f2 = -f1 + 2(f1 dot nhat)nhat
	*   which is equal to f2 = f1 - 2(f1 dot n)n/nlength^2
	*/
	
	// The normal vector to the point of intersection
	double normVec_V = - pos_W*ellipse_V_axis/ellipse_W_axis;
	double normVec_W =   pos_V*ellipse_W_axis/ellipse_V_axis;

	double normVec_length_squared = normVec_V*normVec_V + normVec_W*normVec_W;

	// Dot product of (vel_V0,vel_W0) and the normal vector
	double Vel_dot_NV = vel_V*normVec_V+vel_W*normVec_W;

	// Calculate f2
	double vel_V_2 = -vel_V + 2*Vel_dot_NV*normVec_V/normVec_length_squared;
	double vel_W_2 = -vel_W + 2*Vel_dot_NV*normVec_W/normVec_length_squared;

	// Apply the new velocity vector to the particle globally
	*pvel_V=vel_V_2;
	*pvel_W=vel_W_2;

	// Calculate q and the weighting of the neutron package
	// q=f1-f2
	double delta_vel_V = vel_V-vel_V_2;
	double delta_vel_W = vel_W-vel_W_2;
	double q = V2Q*sqrt( delta_vel_V*delta_vel_V+delta_vel_W*delta_vel_W );

	// Calculate the loss of neutrons due to the reflection
	double mirrorPar[] = {R0, Qc, alpha, Mvalue, W};
	double weight = 1.0;
	StdReflecFunc(q, mirrorPar, &weight);

	return weight;
}

/**
	Use the found side of Intersection to call guide_elliptical_ReflectionOnEllipticSurface with
	the parameters of that side.
*/
double guide_elliptical_handleReflection(double x0, double y0, double z0,
				double *vx_p,double *vy_p,double *vz_p,
				struct SGI *sgi,
				struct Intersection *cc)
{

    if(!sgi->enableSegments){
		if(cc->side == RightSide || cc->side == LeftSide)
			return guide_elliptical_ReflectionOnEllipticSurface(x0,z0,vx_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mArr[cc->side],
					sgi->WArr[cc->side]
					);
		if(cc->side == TopSide || cc->side == BottomSide)
			return guide_elliptical_ReflectionOnEllipticSurface(y0,z0,vy_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mArr[cc->side],
					sgi->WArr[cc->side]
					);
    }
	else{
    	int currentSegment;
    	double combinedLength = 0;
    	int i;
    	for(i=0; i < sgi->numberOfSegments; i++){
    		combinedLength = combinedLength + sgi->segLength[i];
    		if(z0 < combinedLength)	{
    			currentSegment = i; break;
    		}
    	}

		if(cc->side == RightSide)
			return guide_elliptical_ReflectionOnEllipticSurface(x0,z0,vx_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mValuesright[currentSegment],
					sgi->WArr[cc->side] );
		if(cc->side == LeftSide)
			return guide_elliptical_ReflectionOnEllipticSurface(x0,z0,vx_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mValuesleft[currentSegment],
					sgi->WArr[cc->side] );
		if(cc->side == TopSide)
			return guide_elliptical_ReflectionOnEllipticSurface(y0,z0,vy_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mValuestop[currentSegment],
					sgi->WArr[cc->side] );
		if(cc->side == BottomSide)
			return guide_elliptical_ReflectionOnEllipticSurface(y0,z0,vy_p,vz_p,
					sgi->ellipseMinorAxis[cc->side],
					sgi->ellipseMajorAxis[cc->side],
					sgi->ellipseMinorOffset[cc->side],
					sgi->ellipseMajorOffset[cc->side],
					sgi->R0Arr[cc->side],
					sgi->QcArr[cc->side],
					sgi->alphaArr[cc->side],
					sgi->mValuesbottom[currentSegment],
					sgi->WArr[cc->side] );
    }
	return 0;
}

///////////////////////////////////////////////////////////////////////////
/////////////// End of functions
/////////////////////////////////////////////////////////////////////////// 
#line 7700 "./FullInstrument_v4.7_frontend.c"

/* Shared user declarations for all components 'Bender'. */
#line 102 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"

#line 7705 "./FullInstrument_v4.7_frontend.c"

/* Shared user declarations for all components 'PSD_monitor'. */
#line 58 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"

#ifndef ARRAYS_H
#define ARRAYS_H
typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);
#endif
#ifndef ARRAYS_C
#define ARRAYS_C
#include <stdlib.h>

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}
void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}
void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
#endif

#line 7787 "./FullInstrument_v4.7_frontend.c"

/* Shared user declarations for all components 'Monochromator_curved'. */
#line 109 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
#ifndef GAUSS
/* Define these arrays only once for all instances. */
/* Values for Gauss quadrature. Taken from Brice Carnahan, H. A. Luther and
James O Wilkes, "Applied numerical methods", Wiley, 1969, page 103.
This reference is available from the Copenhagen UB2 library */
double Gauss_X[] = {-0.987992518020485, -0.937273392400706, -0.848206583410427,
-0.724417731360170, -0.570972172608539, -0.394151347077563,
-0.201194093997435, 0, 0.201194093997435,
0.394151347077563, 0.570972172608539, 0.724417731360170,
0.848206583410427, 0.937273392400706, 0.987992518020485};
double Gauss_W[] = {0.030753241996117, 0.070366047488108, 0.107159220467172,
0.139570677926154, 0.166269205816994, 0.186161000115562,
0.198431485327111, 0.202578241925561, 0.198431485327111,
0.186161000115562, 0.166269205816994, 0.139570677926154,
0.107159220467172, 0.070366047488108, 0.030753241996117};


#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2*(rms)*(rms)))/(sqrt(2*PI)*(rms)))
#endif


#line 7813 "./FullInstrument_v4.7_frontend.c"

/* Instrument parameters. */
MCNUM mcipslitwidth;
MCNUM mcipvirtual_slit_d;
MCNUM mcipSourceE;
MCNUM mcipDeltaSourceE;
MCNUM mcipEI;
MCNUM mcipA3;
MCNUM mcipA4;
MCNUM mcipSAMPLE;
char* mcipfile_name;
MCNUM mcipRV_mono_Bool;
MCNUM mcipRH_mono_Bool;
MCNUM mcipBe_filter;
MCNUM mcipold_mono;
MCNUM mcipsampleHeight;
MCNUM mcipsampleRadius;
MCNUM mcipsamplePosY;

#define mcNUMIPAR 16
int mcnumipar = 16;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "slitwidth", &mcipslitwidth, instr_type_double, "0.03", 
  "virtual_slit_d", &mcipvirtual_slit_d, instr_type_double, "0.03", 
  "SourceE", &mcipSourceE, instr_type_double, "5", 
  "DeltaSourceE", &mcipDeltaSourceE, instr_type_double, "0.3", 
  "EI", &mcipEI, instr_type_double, "5", 
  "A3", &mcipA3, instr_type_double, "0", 
  "A4", &mcipA4, instr_type_double, "-90", 
  "SAMPLE", &mcipSAMPLE, instr_type_double, "3", 
  "file_name", &mcipfile_name, instr_type_string, "EGECE_13WaveLmax1_geometry.dat.txt", 
  "RV_mono_Bool", &mcipRV_mono_Bool, instr_type_double, "0.5", 
  "RH_mono_Bool", &mcipRH_mono_Bool, instr_type_double, "0.5", 
  "Be_filter", &mcipBe_filter, instr_type_double, "1", 
  "old_mono", &mcipold_mono, instr_type_double, "0", 
  "sampleHeight", &mcipsampleHeight, instr_type_double, "0.01", 
  "sampleRadius", &mcipsampleRadius, instr_type_double, "0.005", 
  "samplePosY", &mcipsamplePosY, instr_type_double, "0.00", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  PSI_CAMEA
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_CAMEA coords_set(0,0,0)
#define slitwidth mcipslitwidth
#define virtual_slit_d mcipvirtual_slit_d
#define SourceE mcipSourceE
#define DeltaSourceE mcipDeltaSourceE
#define EI mcipEI
#define A3 mcipA3
#define A4 mcipA4
#define SAMPLE mcipSAMPLE
#define file_name mcipfile_name
#define RV_mono_Bool mcipRV_mono_Bool
#define RH_mono_Bool mcipRH_mono_Bool
#define Be_filter mcipBe_filter
#define old_mono mcipold_mono
#define sampleHeight mcipsampleHeight
#define sampleRadius mcipsampleRadius
#define samplePosY mcipsamplePosY
#line 42 "FullInstrument_v4.7_frontend.instr"
double slitheight=0.2;
double guide_start=1.5411;
double smallaxis_x_factor5=0.001019;
double smallaxis_y_factor5=0.00193;
double startx5=0.034973;
double starty5=0.11782;
double smallaxis_x_factor3=0.0012272;
double smallaxis_y_factor3=0.11055;
double startx3=0.070026;
double starty3=0.14691;
double lengthfrac2=0.35765;
double smallaxis_x_factor1=0.16453;
double smallaxis_y_factor1=0.093264;
double startx1=0.051571;
double starty1=0.14841;
double length1=13.2556;
double sizeX=0.03;
double sizeY=0.05;
double divreq_x=2;
double divreq_y=2;
double sample_dist=0.03;
double mod_x=0.08;
double mod_y=0.135;
double minimalist_factor=1;
double Mod_sample=35.6;
double closest_element=1.5;
double length4=0.13;
double length5=4.6;
double minlength1=6;
double minlength3=1;
double minlength5=1;
double R05=0.99;
double Qc5=0.0217;
double alpha5=6.07;
double m5=3.5;
double W5=0.003;
double max_smallaxis_x5=0.25;
double max_smallaxis_y5=0.25;
double R03=0.99;
double Qc3=0.0217;
double alpha3=6.07;
double m3=3.5;
double W3=0.003;
double max_smallaxis_x3=0.25;
double max_smallaxis_y3=0.25;
double R02=0.99;
double Qc2=0.0217;
double alpha2=6.07;
double m2=6;
double W2=0.003;
double channels2=1;
double R01=0.99;
double Qc1=0.0217;
double alpha1=6.07;
double m1=3.5;
double W1=0.003;
double max_smallaxis_x1=0.25;
double max_smallaxis_y1=0.25;
double WaveMin=1;
double WaveMax=10;
/********************************************************************************************/
double thetaA1;
double thetaA2;
double R_mono;
double lmin;
double lmax;
/********************************************************************************************/
double u;
double Lambda0;
double dLambda;
double var_divreq_x;
double var_divreq_y;
double x_div;
double y_div;
double startxpoint[8][3][3];
double startypoint[8][3][3];
double startXdirec[8][3];
double startYdirec[8][3];
double startXposition[8][3];
double startYposition[8][3];
double dist;
double position;
double focus_e;
double elength;
double width;
double height;
double curve_small_radius2;
double curveXcenter2[3];
double DeltaA2;
double DeltaB2;
double sinrot3;
double cosrot3;
double rot_factor_g2;
double X1[5];
double X2[5];
double Z1[5];
double Z2[5];
double a[5];
double b[5];
double tmp_double;
double dx_circ;
double dy_circ;
double dr_circ;
double D_circ;
double x1_circ;
double x2_circ;
double y1_circ;
double y2_circ;
double sign_dy;
double k_circ;
double x_solution[2];
double y_solution[2];
double smallaxis_x5;
double smallaxis_y5;
double Linx5;
double Liny5;
double Loutx5;
double Louty5;
double tmp_k;
double tmp_L1;
double tmp_L2;
double tmp_c;
double tmp_b;
double tmp_w1;
double tmp_w2;
double tmp_L;
double endx5;
double endy5;
double endPoint5;
double startx4;
double starty4;
double endx4;
double endy4;
double endPoint4;
double smallaxis_x3;
double smallaxis_y3;
double Linx3;
double Liny3;
double Loutx3;
double Louty3;
double endx3;
double endy3;
double length3;
double minlengthneeded;
double tmp_length;
double endPoint3;
double startx2;
double starty2;
double endx2;
double endy2;
double curve_radius2;
double rot2;
double var_divreq_x_protected;
double var_divreq_y_protected;
double length2;
double endPoint2;
double smallaxis_x1;
double smallaxis_y1;
double Linx1;
double Liny1;
double Loutx1;
double Louty1;
double endx1;
double endy1;
double endPoint1;
double rot_factor_g2;
int flag;
int i;
int part;
int n_check;
int los_logic;
int n1;
int n2;
int line;
int los_tmp[5];
int n_check;
int los_logic;
int los_check;
int ii;
int los_logic_single[7][8];
FILE *fp;
int scat = 0;
// Filter parameters 
double OUTFILTER_THICK;

// Analyser material parameters
double ana_mosaic_h;
double ana_mosaic_v;
double ana_q = 1.87325;
double ana_r0 = 0.8; 
double DA;    // d-spacing analyser
//Analyser geometrical parameters 
double ana_l[9];
//double ana_w=0.01;

// Declarations for 'Coarse Collimator' at the PSD detector surface 
int EntrySlit;
int ExitSlit;
int AnaBlade;
double BladeThickness = 0.0001;// detector coll after 2006, from drawings
double WindowSize = 0.025;
double BladeLength = 0.179;// detector coll after 2006, from drawings
double BladeHeight = 0.272;// detector coll after 2006, from drawings
double FirstWindowSizeL;
double FirstWindowSizeR;
double deltaL;
int coarse;

// Analyzer positions
double z_an1[8];
double y_an1[8];
double thet_an1[8];
int Ana1=0;
double z_an2[8];
double y_an2[8];
double thet_an2[8];
int Ana2=0;
double z_an3[8];
double y_an3[8];
double thet_an3[8];
int Ana3=0;

double pixelstart[8];
double pixelstop[8];

double z_an4[8];
double y_an4[8];
double thet_an4[8];
int Ana4=0;

double z_an5[8];
double y_an5[8];
double thet_an5[8];
int Ana5=0;

double z_an6[8];
double y_an6[8];
double thet_an6[8];
int Ana6=0;

double z_an7[8];
double y_an7[8];
double thet_an7[8];
int Ana7=0;

double z_an8[8];
double y_an8[8];
double thet_an8[8];
int Ana8=0;

// Detector parameters 
double ang_1[8];
double ang_2[8];
double detx_1[8];
double detz_1[8];
double detx_2[8];
double detz_2[8];
double dety_1;
double dety_2;
double det_cen;
int i;

double A3_val;
double A4_val;
double sourceDivergence;

double ANA_MOS_H=60;
double ANA_MOS_V=60;
double AW1=0.012;
double AW2=0.012;
double AW3=0.0125;
double AW4=0.013;
double AW5=0.0135;
double AW6=0.014;
double AW7=0.015;
double AW8=0.016;
double RADCOL_height=0.2;
double RADCOL_length=0.3;
double RADCOL_chan=20;
double RADCOL_radius=0.55;
double RADCOL_foilthick=0.0003;
double RADCOL_thetamin=-5;
double RADCOL_thetamax=5;
int RADCOL_flag=1;
int OUTFILTER=1;

double DET_cen=1.2;
double DET_rad=0.0063;
double DET_len=1.0;
double ana_l_add=0;

int REP=1;
int VIRTUALOUT=0;
int VIRTUALIN=0;
int verbose=0;


double RV_mono;
double RH_mono;

#line 8177 "./FullInstrument_v4.7_frontend.c"
#undef samplePosY
#undef sampleRadius
#undef sampleHeight
#undef old_mono
#undef Be_filter
#undef RH_mono_Bool
#undef RV_mono_Bool
#undef file_name
#undef SAMPLE
#undef A4
#undef A3
#undef EI
#undef DeltaSourceE
#undef SourceE
#undef virtual_slit_d
#undef slitwidth
#undef mcposaPSI_CAMEA
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*23];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[23];
Coords mccomp_posr[23];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[23];
MCNUM  mcPCounter[23];
MCNUM  mcP2Counter[23];
#define mcNUMCOMP 22 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[23];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'Origin' [1]. */
char mccOrigin_profile[16384];
MCNUM mccOrigin_percent;
MCNUM mccOrigin_flag_save;
MCNUM mccOrigin_minutes;

/* Definition parameters for component 'source' [2]. */
#define mccsource_flux_file 0 /* declared as a string. May produce warnings at compile */
#define mccsource_xdiv_file 0 /* declared as a string. May produce warnings at compile */
#define mccsource_ydiv_file 0 /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'source' [2]. */
MCNUM mccsource_radius;
MCNUM mccsource_dist;
MCNUM mccsource_xw;
MCNUM mccsource_yh;
MCNUM mccsource_E0;
MCNUM mccsource_dE;
MCNUM mccsource_Lambda0;
MCNUM mccsource_dLambda;
MCNUM mccsource_I1;
MCNUM mccsource_h;
MCNUM mccsource_w;
MCNUM mccsource_gaussian;
MCNUM mccsource_verbose;
MCNUM mccsource_T1;
MCNUM mccsource_flux_file_perAA;
MCNUM mccsource_flux_file_log;
MCNUM mccsource_Lmin;
MCNUM mccsource_Lmax;
MCNUM mccsource_Emin;
MCNUM mccsource_Emax;
MCNUM mccsource_T2;
MCNUM mccsource_I2;
MCNUM mccsource_T3;
MCNUM mccsource_I3;
MCNUM mccsource_length;
MCNUM mccsource_phi_init;
MCNUM mccsource_theta_init;
MCNUM mccsource_HEtailA;
MCNUM mccsource_HEtailL0;

/* Definition parameters for component 'elliptical_guide_gravity5' [4]. */
#define mccelliptical_guide_gravity5_mvaluesright NULL
#define mccelliptical_guide_gravity5_mvaluesleft NULL
#define mccelliptical_guide_gravity5_mvaluestop NULL
#define mccelliptical_guide_gravity5_mvaluesbottom NULL
#define mccelliptical_guide_gravity5_seglength NULL
/* Setting parameters for component 'elliptical_guide_gravity5' [4]. */
MCNUM mccelliptical_guide_gravity5_l;
MCNUM mccelliptical_guide_gravity5_xwidth;
MCNUM mccelliptical_guide_gravity5_yheight;
MCNUM mccelliptical_guide_gravity5_linxw;
MCNUM mccelliptical_guide_gravity5_loutxw;
MCNUM mccelliptical_guide_gravity5_linyh;
MCNUM mccelliptical_guide_gravity5_loutyh;
MCNUM mccelliptical_guide_gravity5_majorAxisxw;
MCNUM mccelliptical_guide_gravity5_minorAxisxw;
MCNUM mccelliptical_guide_gravity5_majorAxisyh;
MCNUM mccelliptical_guide_gravity5_minorAxisyh;
MCNUM mccelliptical_guide_gravity5_majorAxisoffsetxw;
MCNUM mccelliptical_guide_gravity5_majorAxisoffsetyh;
char mccelliptical_guide_gravity5_dimensionsAt[16384];
char mccelliptical_guide_gravity5_option[16384];
MCNUM mccelliptical_guide_gravity5_R0;
MCNUM mccelliptical_guide_gravity5_Qc;
MCNUM mccelliptical_guide_gravity5_alpha;
MCNUM mccelliptical_guide_gravity5_m;
MCNUM mccelliptical_guide_gravity5_W;
MCNUM mccelliptical_guide_gravity5_alpharight;
MCNUM mccelliptical_guide_gravity5_mright;
MCNUM mccelliptical_guide_gravity5_alphaleft;
MCNUM mccelliptical_guide_gravity5_mleft;
MCNUM mccelliptical_guide_gravity5_alphatop;
MCNUM mccelliptical_guide_gravity5_mtop;
MCNUM mccelliptical_guide_gravity5_alphabottom;
MCNUM mccelliptical_guide_gravity5_mbottom;
char mccelliptical_guide_gravity5_verbose[16384];
MCNUM mccelliptical_guide_gravity5_enableGravity;
MCNUM mccelliptical_guide_gravity5_curvature;

/* Definition parameters for component 'elliptical_guide_gravity3' [7]. */
#define mccelliptical_guide_gravity3_mvaluesright NULL
#define mccelliptical_guide_gravity3_mvaluesleft NULL
#define mccelliptical_guide_gravity3_mvaluestop NULL
#define mccelliptical_guide_gravity3_mvaluesbottom NULL
#define mccelliptical_guide_gravity3_seglength NULL
/* Setting parameters for component 'elliptical_guide_gravity3' [7]. */
MCNUM mccelliptical_guide_gravity3_l;
MCNUM mccelliptical_guide_gravity3_xwidth;
MCNUM mccelliptical_guide_gravity3_yheight;
MCNUM mccelliptical_guide_gravity3_linxw;
MCNUM mccelliptical_guide_gravity3_loutxw;
MCNUM mccelliptical_guide_gravity3_linyh;
MCNUM mccelliptical_guide_gravity3_loutyh;
MCNUM mccelliptical_guide_gravity3_majorAxisxw;
MCNUM mccelliptical_guide_gravity3_minorAxisxw;
MCNUM mccelliptical_guide_gravity3_majorAxisyh;
MCNUM mccelliptical_guide_gravity3_minorAxisyh;
MCNUM mccelliptical_guide_gravity3_majorAxisoffsetxw;
MCNUM mccelliptical_guide_gravity3_majorAxisoffsetyh;
char mccelliptical_guide_gravity3_dimensionsAt[16384];
char mccelliptical_guide_gravity3_option[16384];
MCNUM mccelliptical_guide_gravity3_R0;
MCNUM mccelliptical_guide_gravity3_Qc;
MCNUM mccelliptical_guide_gravity3_alpha;
MCNUM mccelliptical_guide_gravity3_m;
MCNUM mccelliptical_guide_gravity3_W;
MCNUM mccelliptical_guide_gravity3_alpharight;
MCNUM mccelliptical_guide_gravity3_mright;
MCNUM mccelliptical_guide_gravity3_alphaleft;
MCNUM mccelliptical_guide_gravity3_mleft;
MCNUM mccelliptical_guide_gravity3_alphatop;
MCNUM mccelliptical_guide_gravity3_mtop;
MCNUM mccelliptical_guide_gravity3_alphabottom;
MCNUM mccelliptical_guide_gravity3_mbottom;
char mccelliptical_guide_gravity3_verbose[16384];
MCNUM mccelliptical_guide_gravity3_enableGravity;
MCNUM mccelliptical_guide_gravity3_curvature;

/* Setting parameters for component 'Curved_guide_2' [9]. */
MCNUM mccCurved_guide_2_w;
MCNUM mccCurved_guide_2_h;
MCNUM mccCurved_guide_2_r;
MCNUM mccCurved_guide_2_Win;
MCNUM mccCurved_guide_2_k;
MCNUM mccCurved_guide_2_d;
MCNUM mccCurved_guide_2_l;
MCNUM mccCurved_guide_2_R0a;
MCNUM mccCurved_guide_2_Qca;
MCNUM mccCurved_guide_2_alphaa;
MCNUM mccCurved_guide_2_ma;
MCNUM mccCurved_guide_2_Wa;
MCNUM mccCurved_guide_2_R0i;
MCNUM mccCurved_guide_2_Qci;
MCNUM mccCurved_guide_2_alphai;
MCNUM mccCurved_guide_2_mi;
MCNUM mccCurved_guide_2_Wi;
MCNUM mccCurved_guide_2_R0s;
MCNUM mccCurved_guide_2_Qcs;
MCNUM mccCurved_guide_2_alphas;
MCNUM mccCurved_guide_2_ms;
MCNUM mccCurved_guide_2_Ws;

/* Definition parameters for component 'elliptical_guide_gravity1' [11]. */
#define mccelliptical_guide_gravity1_mvaluesright NULL
#define mccelliptical_guide_gravity1_mvaluesleft NULL
#define mccelliptical_guide_gravity1_mvaluestop NULL
#define mccelliptical_guide_gravity1_mvaluesbottom NULL
#define mccelliptical_guide_gravity1_seglength NULL
/* Setting parameters for component 'elliptical_guide_gravity1' [11]. */
MCNUM mccelliptical_guide_gravity1_l;
MCNUM mccelliptical_guide_gravity1_xwidth;
MCNUM mccelliptical_guide_gravity1_yheight;
MCNUM mccelliptical_guide_gravity1_linxw;
MCNUM mccelliptical_guide_gravity1_loutxw;
MCNUM mccelliptical_guide_gravity1_linyh;
MCNUM mccelliptical_guide_gravity1_loutyh;
MCNUM mccelliptical_guide_gravity1_majorAxisxw;
MCNUM mccelliptical_guide_gravity1_minorAxisxw;
MCNUM mccelliptical_guide_gravity1_majorAxisyh;
MCNUM mccelliptical_guide_gravity1_minorAxisyh;
MCNUM mccelliptical_guide_gravity1_majorAxisoffsetxw;
MCNUM mccelliptical_guide_gravity1_majorAxisoffsetyh;
char mccelliptical_guide_gravity1_dimensionsAt[16384];
char mccelliptical_guide_gravity1_option[16384];
MCNUM mccelliptical_guide_gravity1_R0;
MCNUM mccelliptical_guide_gravity1_Qc;
MCNUM mccelliptical_guide_gravity1_alpha;
MCNUM mccelliptical_guide_gravity1_m;
MCNUM mccelliptical_guide_gravity1_W;
MCNUM mccelliptical_guide_gravity1_alpharight;
MCNUM mccelliptical_guide_gravity1_mright;
MCNUM mccelliptical_guide_gravity1_alphaleft;
MCNUM mccelliptical_guide_gravity1_mleft;
MCNUM mccelliptical_guide_gravity1_alphatop;
MCNUM mccelliptical_guide_gravity1_mtop;
MCNUM mccelliptical_guide_gravity1_alphabottom;
MCNUM mccelliptical_guide_gravity1_mbottom;
char mccelliptical_guide_gravity1_verbose[16384];
MCNUM mccelliptical_guide_gravity1_enableGravity;
MCNUM mccelliptical_guide_gravity1_curvature;

/* Setting parameters for component 'SlitGuideEnd' [13]. */
MCNUM mccSlitGuideEnd_xmin;
MCNUM mccSlitGuideEnd_xmax;
MCNUM mccSlitGuideEnd_ymin;
MCNUM mccSlitGuideEnd_ymax;
MCNUM mccSlitGuideEnd_radius;
MCNUM mccSlitGuideEnd_xwidth;
MCNUM mccSlitGuideEnd_yheight;

/* Setting parameters for component 'SlitMonitor' [14]. */
int mccSlitMonitor_nx;
int mccSlitMonitor_ny;
char mccSlitMonitor_filename[16384];
MCNUM mccSlitMonitor_xmin;
MCNUM mccSlitMonitor_xmax;
MCNUM mccSlitMonitor_ymin;
MCNUM mccSlitMonitor_ymax;
MCNUM mccSlitMonitor_xwidth;
MCNUM mccSlitMonitor_yheight;
MCNUM mccSlitMonitor_restore_neutron;
int mccSlitMonitor_nowritefile;

/* Setting parameters for component 'MonoFLat' [15]. */
char mccMonoFLat_reflect[16384];
char mccMonoFLat_transmit[16384];
MCNUM mccMonoFLat_zwidth;
MCNUM mccMonoFLat_yheight;
MCNUM mccMonoFLat_gap;
MCNUM mccMonoFLat_NH;
MCNUM mccMonoFLat_NV;
MCNUM mccMonoFLat_mosaich;
MCNUM mccMonoFLat_mosaicv;
MCNUM mccMonoFLat_r0;
MCNUM mccMonoFLat_t0;
MCNUM mccMonoFLat_Q;
MCNUM mccMonoFLat_RV;
MCNUM mccMonoFLat_RH;
MCNUM mccMonoFLat_DM;
MCNUM mccMonoFLat_mosaic;
MCNUM mccMonoFLat_width;
MCNUM mccMonoFLat_height;
MCNUM mccMonoFLat_verbose;
MCNUM mccMonoFLat_order;

/* Setting parameters for component 'Mono' [16]. */
char mccMono_reflect[16384];
char mccMono_transmit[16384];
MCNUM mccMono_zwidth;
MCNUM mccMono_yheight;
MCNUM mccMono_gap;
MCNUM mccMono_NH;
MCNUM mccMono_NV;
MCNUM mccMono_mosaich;
MCNUM mccMono_mosaicv;
MCNUM mccMono_r0;
MCNUM mccMono_t0;
MCNUM mccMono_Q;
MCNUM mccMono_RV;
MCNUM mccMono_RH;
MCNUM mccMono_DM;
MCNUM mccMono_mosaic;
MCNUM mccMono_width;
MCNUM mccMono_height;
MCNUM mccMono_verbose;
MCNUM mccMono_order;

/* Setting parameters for component 'slit' [18]. */
MCNUM mccslit_xmin;
MCNUM mccslit_xmax;
MCNUM mccslit_ymin;
MCNUM mccslit_ymax;
MCNUM mccslit_radius;
MCNUM mccslit_xwidth;
MCNUM mccslit_yheight;

/* Setting parameters for component 'psd_monitor_sample_pos' [20]. */
int mccpsd_monitor_sample_pos_nx;
int mccpsd_monitor_sample_pos_ny;
char mccpsd_monitor_sample_pos_filename[16384];
MCNUM mccpsd_monitor_sample_pos_xmin;
MCNUM mccpsd_monitor_sample_pos_xmax;
MCNUM mccpsd_monitor_sample_pos_ymin;
MCNUM mccpsd_monitor_sample_pos_ymax;
MCNUM mccpsd_monitor_sample_pos_xwidth;
MCNUM mccpsd_monitor_sample_pos_yheight;
MCNUM mccpsd_monitor_sample_pos_restore_neutron;
int mccpsd_monitor_sample_pos_nowritefile;

/* Definition parameters for component 'e_monitor_sample_ps' [21]. */
#define mcce_monitor_sample_ps_nE 200
/* Setting parameters for component 'e_monitor_sample_ps' [21]. */
char mcce_monitor_sample_ps_filename[16384];
MCNUM mcce_monitor_sample_ps_xmin;
MCNUM mcce_monitor_sample_ps_xmax;
MCNUM mcce_monitor_sample_ps_ymin;
MCNUM mcce_monitor_sample_ps_ymax;
MCNUM mcce_monitor_sample_ps_xwidth;
MCNUM mcce_monitor_sample_ps_yheight;
MCNUM mcce_monitor_sample_ps_Emin;
MCNUM mcce_monitor_sample_ps_Emax;
MCNUM mcce_monitor_sample_ps_restore_neutron;
int mcce_monitor_sample_ps_nowritefile;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 44 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 8528 "./FullInstrument_v4.7_frontend.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source' [2]. */
#define mccompcurname  source
#define mccompcurtype  Source_gen4
#define mccompcurindex 2
#define flux_file mccsource_flux_file
#define xdiv_file mccsource_xdiv_file
#define ydiv_file mccsource_ydiv_file
#define p_in mccsource_p_in
#define lambda0 mccsource_lambda0
#define lambda02 mccsource_lambda02
#define L2P mccsource_L2P
#define lambda0b mccsource_lambda0b
#define lambda02b mccsource_lambda02b
#define L2Pb mccsource_L2Pb
#define lambda0c mccsource_lambda0c
#define lambda02c mccsource_lambda02c
#define L2Pc mccsource_L2Pc
#define pTable mccsource_pTable
#define pTable_x mccsource_pTable_x
#define pTable_y mccsource_pTable_y
#define radius mccsource_radius
#define dist mccsource_dist
#define xw mccsource_xw
#define yh mccsource_yh
#define E0 mccsource_E0
#define dE mccsource_dE
#define Lambda0 mccsource_Lambda0
#define dLambda mccsource_dLambda
#define I1 mccsource_I1
#define h mccsource_h
#define w mccsource_w
#define gaussian mccsource_gaussian
#define verbose mccsource_verbose
#define T1 mccsource_T1
#define flux_file_perAA mccsource_flux_file_perAA
#define flux_file_log mccsource_flux_file_log
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define Emin mccsource_Emin
#define Emax mccsource_Emax
#define T2 mccsource_T2
#define I2 mccsource_I2
#define T3 mccsource_T3
#define I3 mccsource_I3
#define length mccsource_length
#define phi_init mccsource_phi_init
#define theta_init mccsource_theta_init
#define HEtailA mccsource_HEtailA
#define HEtailL0 mccsource_HEtailL0
#line 169 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
  double p_in;
  double Lambda0;
  double lambda0, lambda02, L2P;     /* first Maxwellian source */
  double lambda0b, lambda02b, L2Pb;  /* second Maxwellian source */
  double lambda0c, lambda02c, L2Pc;  /* third Maxwellian source */
  t_Table pTable;
  t_Table pTable_x;
  t_Table pTable_y;
  double  pTable_xmin;
  double  pTable_xmax;
  double  pTable_xsum=0;
  double  pTable_ymin;
  double  pTable_ymax;
  double  pTable_ysum=0;
  double  pTable_dxmin;
  double  pTable_dxmax;
  double  pTable_dymin;
  double  pTable_dymax;
#line 8609 "./FullInstrument_v4.7_frontend.c"
#undef HEtailL0
#undef HEtailA
#undef theta_init
#undef phi_init
#undef length
#undef I3
#undef T3
#undef I2
#undef T2
#undef Emax
#undef Emin
#undef Lmax
#undef Lmin
#undef flux_file_log
#undef flux_file_perAA
#undef T1
#undef verbose
#undef gaussian
#undef w
#undef h
#undef I1
#undef dLambda
#undef Lambda0
#undef dE
#undef E0
#undef yh
#undef xw
#undef dist
#undef radius
#undef pTable_y
#undef pTable_x
#undef pTable
#undef L2Pc
#undef lambda02c
#undef lambda0c
#undef L2Pb
#undef lambda02b
#undef lambda0b
#undef L2P
#undef lambda02
#undef lambda0
#undef p_in
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'StartOfGuide' [3]. */
#define mccompcurname  StartOfGuide
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'elliptical_guide_gravity5' [4]. */
#define mccompcurname  elliptical_guide_gravity5
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 4
#define mvaluesright mccelliptical_guide_gravity5_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity5_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity5_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity5_mvaluesbottom
#define seglength mccelliptical_guide_gravity5_seglength
#define guideInfo mccelliptical_guide_gravity5_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity5_latestParticleCollision
#define Gx mccelliptical_guide_gravity5_Gx
#define Gy mccelliptical_guide_gravity5_Gy
#define Gz mccelliptical_guide_gravity5_Gz
#define Gx0 mccelliptical_guide_gravity5_Gx0
#define Gy0 mccelliptical_guide_gravity5_Gy0
#define Gz0 mccelliptical_guide_gravity5_Gz0
#define Circ mccelliptical_guide_gravity5_Circ
#define dynamicalSegLength mccelliptical_guide_gravity5_dynamicalSegLength
#define l mccelliptical_guide_gravity5_l
#define xwidth mccelliptical_guide_gravity5_xwidth
#define yheight mccelliptical_guide_gravity5_yheight
#define linxw mccelliptical_guide_gravity5_linxw
#define loutxw mccelliptical_guide_gravity5_loutxw
#define linyh mccelliptical_guide_gravity5_linyh
#define loutyh mccelliptical_guide_gravity5_loutyh
#define majorAxisxw mccelliptical_guide_gravity5_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity5_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity5_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity5_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity5_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity5_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity5_dimensionsAt
#define option mccelliptical_guide_gravity5_option
#define R0 mccelliptical_guide_gravity5_R0
#define Qc mccelliptical_guide_gravity5_Qc
#define alpha mccelliptical_guide_gravity5_alpha
#define m mccelliptical_guide_gravity5_m
#define W mccelliptical_guide_gravity5_W
#define alpharight mccelliptical_guide_gravity5_alpharight
#define mright mccelliptical_guide_gravity5_mright
#define alphaleft mccelliptical_guide_gravity5_alphaleft
#define mleft mccelliptical_guide_gravity5_mleft
#define alphatop mccelliptical_guide_gravity5_alphatop
#define mtop mccelliptical_guide_gravity5_mtop
#define alphabottom mccelliptical_guide_gravity5_alphabottom
#define mbottom mccelliptical_guide_gravity5_mbottom
#define verbose mccelliptical_guide_gravity5_verbose
#define enableGravity mccelliptical_guide_gravity5_enableGravity
#define curvature mccelliptical_guide_gravity5_curvature
#line 803 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
	/** 
		All variables below is locally declared 
		and hence accessible through OUTPUT PARAMETERS.
	*/ 
	struct SGI guideInfo; // Static Guide information, is set in INITIALIZE
	struct Intersection latestParticleCollision; // Is changed duing trace
	double Gx,Gy,Gz;  // Local gravity vector, is set once in INITIALIZE
	double Gx0, Gy0, Gz0;
	double Circ;
	double *dynamicalSegLength;

#line 8729 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EndOfelement_5' [5]. */
#define mccompcurname  EndOfelement_5
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EndOfelement_4' [6]. */
#define mccompcurname  EndOfelement_4
#define mccompcurtype  Arm
#define mccompcurindex 6
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'elliptical_guide_gravity3' [7]. */
#define mccompcurname  elliptical_guide_gravity3
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 7
#define mvaluesright mccelliptical_guide_gravity3_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity3_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity3_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity3_mvaluesbottom
#define seglength mccelliptical_guide_gravity3_seglength
#define guideInfo mccelliptical_guide_gravity3_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity3_latestParticleCollision
#define Gx mccelliptical_guide_gravity3_Gx
#define Gy mccelliptical_guide_gravity3_Gy
#define Gz mccelliptical_guide_gravity3_Gz
#define Gx0 mccelliptical_guide_gravity3_Gx0
#define Gy0 mccelliptical_guide_gravity3_Gy0
#define Gz0 mccelliptical_guide_gravity3_Gz0
#define Circ mccelliptical_guide_gravity3_Circ
#define dynamicalSegLength mccelliptical_guide_gravity3_dynamicalSegLength
#define l mccelliptical_guide_gravity3_l
#define xwidth mccelliptical_guide_gravity3_xwidth
#define yheight mccelliptical_guide_gravity3_yheight
#define linxw mccelliptical_guide_gravity3_linxw
#define loutxw mccelliptical_guide_gravity3_loutxw
#define linyh mccelliptical_guide_gravity3_linyh
#define loutyh mccelliptical_guide_gravity3_loutyh
#define majorAxisxw mccelliptical_guide_gravity3_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity3_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity3_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity3_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity3_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity3_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity3_dimensionsAt
#define option mccelliptical_guide_gravity3_option
#define R0 mccelliptical_guide_gravity3_R0
#define Qc mccelliptical_guide_gravity3_Qc
#define alpha mccelliptical_guide_gravity3_alpha
#define m mccelliptical_guide_gravity3_m
#define W mccelliptical_guide_gravity3_W
#define alpharight mccelliptical_guide_gravity3_alpharight
#define mright mccelliptical_guide_gravity3_mright
#define alphaleft mccelliptical_guide_gravity3_alphaleft
#define mleft mccelliptical_guide_gravity3_mleft
#define alphatop mccelliptical_guide_gravity3_alphatop
#define mtop mccelliptical_guide_gravity3_mtop
#define alphabottom mccelliptical_guide_gravity3_alphabottom
#define mbottom mccelliptical_guide_gravity3_mbottom
#define verbose mccelliptical_guide_gravity3_verbose
#define enableGravity mccelliptical_guide_gravity3_enableGravity
#define curvature mccelliptical_guide_gravity3_curvature
#line 803 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
	/** 
		All variables below is locally declared 
		and hence accessible through OUTPUT PARAMETERS.
	*/ 
	struct SGI guideInfo; // Static Guide information, is set in INITIALIZE
	struct Intersection latestParticleCollision; // Is changed duing trace
	double Gx,Gy,Gz;  // Local gravity vector, is set once in INITIALIZE
	double Gx0, Gy0, Gz0;
	double Circ;
	double *dynamicalSegLength;

#line 8858 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EndOfelement_3' [8]. */
#define mccompcurname  EndOfelement_3
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Curved_guide_2' [9]. */
#define mccompcurname  Curved_guide_2
#define mccompcurtype  Bender
#define mccompcurindex 9
#define bk mccCurved_guide_2_bk
#define mWin mccCurved_guide_2_mWin
#define w mccCurved_guide_2_w
#define h mccCurved_guide_2_h
#define r mccCurved_guide_2_r
#define Win mccCurved_guide_2_Win
#define k mccCurved_guide_2_k
#define d mccCurved_guide_2_d
#define l mccCurved_guide_2_l
#define R0a mccCurved_guide_2_R0a
#define Qca mccCurved_guide_2_Qca
#define alphaa mccCurved_guide_2_alphaa
#define ma mccCurved_guide_2_ma
#define Wa mccCurved_guide_2_Wa
#define R0i mccCurved_guide_2_R0i
#define Qci mccCurved_guide_2_Qci
#define alphai mccCurved_guide_2_alphai
#define mi mccCurved_guide_2_mi
#define Wi mccCurved_guide_2_Wi
#define R0s mccCurved_guide_2_R0s
#define Qcs mccCurved_guide_2_Qcs
#define alphas mccCurved_guide_2_alphas
#define ms mccCurved_guide_2_ms
#define Ws mccCurved_guide_2_Ws
#line 108 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
double bk, mWin;
#line 8947 "./FullInstrument_v4.7_frontend.c"
#undef Ws
#undef ms
#undef alphas
#undef Qcs
#undef R0s
#undef Wi
#undef mi
#undef alphai
#undef Qci
#undef R0i
#undef Wa
#undef ma
#undef alphaa
#undef Qca
#undef R0a
#undef l
#undef d
#undef k
#undef Win
#undef r
#undef h
#undef w
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EndOfelement_2' [10]. */
#define mccompcurname  EndOfelement_2
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'elliptical_guide_gravity1' [11]. */
#define mccompcurname  elliptical_guide_gravity1
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 11
#define mvaluesright mccelliptical_guide_gravity1_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity1_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity1_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity1_mvaluesbottom
#define seglength mccelliptical_guide_gravity1_seglength
#define guideInfo mccelliptical_guide_gravity1_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity1_latestParticleCollision
#define Gx mccelliptical_guide_gravity1_Gx
#define Gy mccelliptical_guide_gravity1_Gy
#define Gz mccelliptical_guide_gravity1_Gz
#define Gx0 mccelliptical_guide_gravity1_Gx0
#define Gy0 mccelliptical_guide_gravity1_Gy0
#define Gz0 mccelliptical_guide_gravity1_Gz0
#define Circ mccelliptical_guide_gravity1_Circ
#define dynamicalSegLength mccelliptical_guide_gravity1_dynamicalSegLength
#define l mccelliptical_guide_gravity1_l
#define xwidth mccelliptical_guide_gravity1_xwidth
#define yheight mccelliptical_guide_gravity1_yheight
#define linxw mccelliptical_guide_gravity1_linxw
#define loutxw mccelliptical_guide_gravity1_loutxw
#define linyh mccelliptical_guide_gravity1_linyh
#define loutyh mccelliptical_guide_gravity1_loutyh
#define majorAxisxw mccelliptical_guide_gravity1_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity1_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity1_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity1_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity1_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity1_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity1_dimensionsAt
#define option mccelliptical_guide_gravity1_option
#define R0 mccelliptical_guide_gravity1_R0
#define Qc mccelliptical_guide_gravity1_Qc
#define alpha mccelliptical_guide_gravity1_alpha
#define m mccelliptical_guide_gravity1_m
#define W mccelliptical_guide_gravity1_W
#define alpharight mccelliptical_guide_gravity1_alpharight
#define mright mccelliptical_guide_gravity1_mright
#define alphaleft mccelliptical_guide_gravity1_alphaleft
#define mleft mccelliptical_guide_gravity1_mleft
#define alphatop mccelliptical_guide_gravity1_alphatop
#define mtop mccelliptical_guide_gravity1_mtop
#define alphabottom mccelliptical_guide_gravity1_alphabottom
#define mbottom mccelliptical_guide_gravity1_mbottom
#define verbose mccelliptical_guide_gravity1_verbose
#define enableGravity mccelliptical_guide_gravity1_enableGravity
#define curvature mccelliptical_guide_gravity1_curvature
#line 803 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
	/** 
		All variables below is locally declared 
		and hence accessible through OUTPUT PARAMETERS.
	*/ 
	struct SGI guideInfo; // Static Guide information, is set in INITIALIZE
	struct Intersection latestParticleCollision; // Is changed duing trace
	double Gx,Gy,Gz;  // Local gravity vector, is set once in INITIALIZE
	double Gx0, Gy0, Gz0;
	double Circ;
	double *dynamicalSegLength;

#line 9046 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'EndOfelement_1' [12]. */
#define mccompcurname  EndOfelement_1
#define mccompcurtype  Arm
#define mccompcurindex 12
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'SlitGuideEnd' [13]. */
#define mccompcurname  SlitGuideEnd
#define mccompcurtype  Slit
#define mccompcurindex 13
#define xmin mccSlitGuideEnd_xmin
#define xmax mccSlitGuideEnd_xmax
#define ymin mccSlitGuideEnd_ymin
#define ymax mccSlitGuideEnd_ymax
#define radius mccSlitGuideEnd_radius
#define xwidth mccSlitGuideEnd_xwidth
#define yheight mccSlitGuideEnd_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'SlitMonitor' [14]. */
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
#define nx mccSlitMonitor_nx
#define ny mccSlitMonitor_ny
#define filename mccSlitMonitor_filename
#define xmin mccSlitMonitor_xmin
#define xmax mccSlitMonitor_xmax
#define ymin mccSlitMonitor_ymin
#define ymax mccSlitMonitor_ymax
#define xwidth mccSlitMonitor_xwidth
#define yheight mccSlitMonitor_yheight
#define restore_neutron mccSlitMonitor_restore_neutron
#define nowritefile mccSlitMonitor_nowritefile
#line 63 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9149 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'MonoFLat' [15]. */
#define mccompcurname  MonoFLat
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 15
#define mos_rms_y mccMonoFLat_mos_rms_y
#define mos_rms_z mccMonoFLat_mos_rms_z
#define mos_rms_max mccMonoFLat_mos_rms_max
#define mono_Q mccMonoFLat_mono_Q
#define SlabWidth mccMonoFLat_SlabWidth
#define SlabHeight mccMonoFLat_SlabHeight
#define rTable mccMonoFLat_rTable
#define tTable mccMonoFLat_tTable
#define row mccMonoFLat_row
#define col mccMonoFLat_col
#define tiltH mccMonoFLat_tiltH
#define tiltV mccMonoFLat_tiltV
#define reflect mccMonoFLat_reflect
#define transmit mccMonoFLat_transmit
#define zwidth mccMonoFLat_zwidth
#define yheight mccMonoFLat_yheight
#define gap mccMonoFLat_gap
#define NH mccMonoFLat_NH
#define NV mccMonoFLat_NV
#define mosaich mccMonoFLat_mosaich
#define mosaicv mccMonoFLat_mosaicv
#define r0 mccMonoFLat_r0
#define t0 mccMonoFLat_t0
#define Q mccMonoFLat_Q
#define RV mccMonoFLat_RV
#define RH mccMonoFLat_RH
#define DM mccMonoFLat_DM
#define mosaic mccMonoFLat_mosaic
#define width mccMonoFLat_width
#define height mccMonoFLat_height
#define verbose mccMonoFLat_verbose
#define order mccMonoFLat_order
#line 136 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
  double SlabWidth, SlabHeight;
  t_Table rTable, tTable;
  double row,col;
  double* tiltH;
  double* tiltV;
#line 9214 "./FullInstrument_v4.7_frontend.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Mono' [16]. */
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 16
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define SlabWidth mccMono_SlabWidth
#define SlabHeight mccMono_SlabHeight
#define rTable mccMono_rTable
#define tTable mccMono_tTable
#define row mccMono_row
#define col mccMono_col
#define tiltH mccMono_tiltH
#define tiltV mccMono_tiltV
#define reflect mccMono_reflect
#define transmit mccMono_transmit
#define zwidth mccMono_zwidth
#define yheight mccMono_yheight
#define gap mccMono_gap
#define NH mccMono_NH
#define NV mccMono_NV
#define mosaich mccMono_mosaich
#define mosaicv mccMono_mosaicv
#define r0 mccMono_r0
#define t0 mccMono_t0
#define Q mccMono_Q
#define RV mccMono_RV
#define RH mccMono_RH
#define DM mccMono_DM
#define mosaic mccMono_mosaic
#define width mccMono_width
#define height mccMono_height
#define verbose mccMono_verbose
#define order mccMono_order
#line 136 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
  double mos_rms_y; /* root-mean-square of mosaic, in radians */
  double mos_rms_z;
  double mos_rms_max;
  double mono_Q;
  double SlabWidth, SlabHeight;
  t_Table rTable, tTable;
  double row,col;
  double* tiltH;
  double* tiltV;
#line 9297 "./FullInstrument_v4.7_frontend.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'A2' [17]. */
#define mccompcurname  A2
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slit' [18]. */
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccslit_xmin
#define xmax mccslit_xmax
#define ymin mccslit_ymin
#define ymax mccslit_ymax
#define radius mccslit_radius
#define xwidth mccslit_xwidth
#define yheight mccslit_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'SamplePos' [19]. */
#define mccompcurname  SamplePos
#define mccompcurtype  Arm
#define mccompcurindex 19
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_sample_pos' [20]. */
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
#define nx mccpsd_monitor_sample_pos_nx
#define ny mccpsd_monitor_sample_pos_ny
#define filename mccpsd_monitor_sample_pos_filename
#define xmin mccpsd_monitor_sample_pos_xmin
#define xmax mccpsd_monitor_sample_pos_xmax
#define ymin mccpsd_monitor_sample_pos_ymin
#define ymax mccpsd_monitor_sample_pos_ymax
#define xwidth mccpsd_monitor_sample_pos_xwidth
#define yheight mccpsd_monitor_sample_pos_yheight
#define restore_neutron mccpsd_monitor_sample_pos_restore_neutron
#define nowritefile mccpsd_monitor_sample_pos_nowritefile
#line 63 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9394 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'e_monitor_sample_ps' [21]. */
#define mccompcurname  e_monitor_sample_ps
#define mccompcurtype  E_monitor
#define mccompcurindex 21
#define nE mcce_monitor_sample_ps_nE
#define E_N mcce_monitor_sample_ps_E_N
#define E_p mcce_monitor_sample_ps_E_p
#define E_p2 mcce_monitor_sample_ps_E_p2
#define S_p mcce_monitor_sample_ps_S_p
#define S_pE mcce_monitor_sample_ps_S_pE
#define S_pE2 mcce_monitor_sample_ps_S_pE2
#define filename mcce_monitor_sample_ps_filename
#define xmin mcce_monitor_sample_ps_xmin
#define xmax mcce_monitor_sample_ps_xmax
#define ymin mcce_monitor_sample_ps_ymin
#define ymax mcce_monitor_sample_ps_ymax
#define xwidth mcce_monitor_sample_ps_xwidth
#define yheight mcce_monitor_sample_ps_yheight
#define Emin mcce_monitor_sample_ps_Emin
#define Emax mcce_monitor_sample_ps_Emax
#define restore_neutron mcce_monitor_sample_ps_restore_neutron
#define nowritefile mcce_monitor_sample_ps_nowritefile
#line 60 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
double E_N[nE];
double E_p[nE], E_p2[nE];
double S_p, S_pE, S_pE2;
#line 9439 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposaStartOfGuide, mcposrStartOfGuide;
Rotation mcrotaStartOfGuide, mcrotrStartOfGuide;
Coords mcposaelliptical_guide_gravity5, mcposrelliptical_guide_gravity5;
Rotation mcrotaelliptical_guide_gravity5, mcrotrelliptical_guide_gravity5;
Coords mcposaEndOfelement_5, mcposrEndOfelement_5;
Rotation mcrotaEndOfelement_5, mcrotrEndOfelement_5;
Coords mcposaEndOfelement_4, mcposrEndOfelement_4;
Rotation mcrotaEndOfelement_4, mcrotrEndOfelement_4;
Coords mcposaelliptical_guide_gravity3, mcposrelliptical_guide_gravity3;
Rotation mcrotaelliptical_guide_gravity3, mcrotrelliptical_guide_gravity3;
Coords mcposaEndOfelement_3, mcposrEndOfelement_3;
Rotation mcrotaEndOfelement_3, mcrotrEndOfelement_3;
Coords mcposaCurved_guide_2, mcposrCurved_guide_2;
Rotation mcrotaCurved_guide_2, mcrotrCurved_guide_2;
Coords mcposaEndOfelement_2, mcposrEndOfelement_2;
Rotation mcrotaEndOfelement_2, mcrotrEndOfelement_2;
Coords mcposaelliptical_guide_gravity1, mcposrelliptical_guide_gravity1;
Rotation mcrotaelliptical_guide_gravity1, mcrotrelliptical_guide_gravity1;
Coords mcposaEndOfelement_1, mcposrEndOfelement_1;
Rotation mcrotaEndOfelement_1, mcrotrEndOfelement_1;
Coords mcposaSlitGuideEnd, mcposrSlitGuideEnd;
Rotation mcrotaSlitGuideEnd, mcrotrSlitGuideEnd;
Coords mcposaSlitMonitor, mcposrSlitMonitor;
Rotation mcrotaSlitMonitor, mcrotrSlitMonitor;
Coords mcposaMonoFLat, mcposrMonoFLat;
Rotation mcrotaMonoFLat, mcrotrMonoFLat;
Coords mcposaMono, mcposrMono;
Rotation mcrotaMono, mcrotrMono;
Coords mcposaA2, mcposrA2;
Rotation mcrotaA2, mcrotrA2;
Coords mcposaslit, mcposrslit;
Rotation mcrotaslit, mcrotrslit;
Coords mcposaSamplePos, mcposrSamplePos;
Rotation mcrotaSamplePos, mcrotrSamplePos;
Coords mcposapsd_monitor_sample_pos, mcposrpsd_monitor_sample_pos;
Rotation mcrotapsd_monitor_sample_pos, mcrotrpsd_monitor_sample_pos;
Coords mcposae_monitor_sample_ps, mcposre_monitor_sample_ps;
Rotation mcrotae_monitor_sample_ps, mcrotre_monitor_sample_ps;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  PSI_CAMEA
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaPSI_CAMEA coords_set(0,0,0)
#define slitwidth mcipslitwidth
#define virtual_slit_d mcipvirtual_slit_d
#define SourceE mcipSourceE
#define DeltaSourceE mcipDeltaSourceE
#define EI mcipEI
#define A3 mcipA3
#define A4 mcipA4
#define SAMPLE mcipSAMPLE
#define file_name mcipfile_name
#define RV_mono_Bool mcipRV_mono_Bool
#define RH_mono_Bool mcipRH_mono_Bool
#define Be_filter mcipBe_filter
#define old_mono mcipold_mono
#define sampleHeight mcipsampleHeight
#define sampleRadius mcipsampleRadius
#define samplePosY mcipsamplePosY
#line 344 "FullInstrument_v4.7_frontend.instr"
{
/********************************************************************************************/
thetaA1= asin(1.0/(0.1106*sqrt(EI)*2.0*3.355));
thetaA2=2.0*thetaA1;
R_mono=sqrt( ( 2.0*1.6*1.6+2.0*1.6*1.6*cos(thetaA2)) )/2.0/sqrt(sin(thetaA2)*sin(thetaA2))*2;
thetaA2=thetaA2*RAD2DEG;
thetaA1=thetaA1*RAD2DEG;
lmin=1.0/0.1106/sqrt(SourceE+DeltaSourceE/2.0);
lmax=1.0/0.1106/sqrt(SourceE-DeltaSourceE/2.0);
printf("\n lmin: %f      lmax: %f \n",lmin,lmax);
printf("\n E0: %f      E0max: %f E0min %f \n",SourceE,SourceE+DeltaSourceE/2.0,SourceE-DeltaSourceE/2.0);
/********************************************************************************************/

dLambda = 0.5*(WaveMax - WaveMin);
Lambda0 = dLambda+WaveMin;

var_divreq_x = divreq_x;
var_divreq_y = divreq_y;

u=1e-4;


endPoint1 = Mod_sample - sample_dist;

endy1 = sizeY + 2*sample_dist*tan(divreq_y*DEG2RAD);


endx1 = sizeX + 2*sample_dist*tan(divreq_x*DEG2RAD);


var_divreq_x = atan(endx1*tan(var_divreq_x*DEG2RAD)/startx1)*RAD2DEG;
var_divreq_y = atan(endy1*tan(var_divreq_y*DEG2RAD)/starty1)*RAD2DEG;



if (startx1 > endx1) {
smallaxis_x1 = 0.5*startx1+ 0.5*smallaxis_x_factor1*(max_smallaxis_x1-startx1);
}
else {
smallaxis_x1 = 0.5*endx1+ 0.5*smallaxis_x_factor1*(max_smallaxis_x1-endx1);
}
if (starty1 > endy1) {
smallaxis_y1 = 0.5*starty1+ 0.5*smallaxis_y_factor1*(max_smallaxis_y1-starty1);
}
else {
smallaxis_y1 = 0.5*endy1+ 0.5*smallaxis_y_factor1*(max_smallaxis_y1-endy1);
}


// calculating focus points for elliptic guide 1 x direction
tmp_w1=startx1;
tmp_w2=endx1;
tmp_b=smallaxis_x1;
tmp_L=length1;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Linx1=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Loutx1= Linx1+tmp_L1-tmp_L2;

// calculating focus points for elliptic guide 1 y direction
tmp_w1=starty1;
tmp_w2=endy1;
tmp_b=smallaxis_y1;
tmp_L=length1;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Liny1=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Louty1= Liny1+tmp_L1-tmp_L2;



endPoint2 = endPoint1 - length1 - u;

minlengthneeded = 0;



length2 = minlengthneeded + (endPoint2 - guide_start - minlengthneeded - length4 - length5 - minlength3 - minlength5)*lengthfrac2;


var_divreq_x_protected = var_divreq_x;
var_divreq_y_protected = var_divreq_y;

rot2 = 0;
los_logic = 1;
while(los_logic==1) {
rot2 = rot2 + 0.0002;
curve_radius2 = length2/(rot2*DEG2RAD);

var_divreq_x = var_divreq_x_protected;
var_divreq_y = var_divreq_y_protected;



endx2 = startx1;
endy2 = starty1;


startx2 = endx2;
starty2 = endy2;


endPoint3 = endPoint2 - length2 - u;

minlengthneeded = 0;

if (minlengthneeded < 1) minlengthneeded = 1;



length3 = endPoint3 - guide_start - length4 - length5;

endx3 = startx2;
endy3 = starty2;


var_divreq_x = atan(endx3*tan(var_divreq_x*DEG2RAD)/startx3)*RAD2DEG;
var_divreq_y = atan(endy3*tan(var_divreq_y*DEG2RAD)/starty3)*RAD2DEG;



if (startx3 > endx3) {
smallaxis_x3 = 0.5*startx3+ 0.5*smallaxis_x_factor3*(max_smallaxis_x3-startx3);
}
else {
smallaxis_x3 = 0.5*endx3+ 0.5*smallaxis_x_factor3*(max_smallaxis_x3-endx3);
}
if (starty3 > endy3) {
smallaxis_y3 = 0.5*starty3+ 0.5*smallaxis_y_factor3*(max_smallaxis_y3-starty3);
}
else {
smallaxis_y3 = 0.5*endy3+ 0.5*smallaxis_y_factor3*(max_smallaxis_y3-endy3);
}


// calculating focus points for elliptic guide 3 x direction
tmp_w1=startx3;
tmp_w2=endx3;
tmp_b=smallaxis_x3;
tmp_L=length3;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Linx3=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Loutx3= Linx3+tmp_L1-tmp_L2;

// calculating focus points for elliptic guide 3 y direction
tmp_w1=starty3;
tmp_w2=endy3;
tmp_b=smallaxis_y3;
tmp_L=length3;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Liny3=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Louty3= Liny3+tmp_L1-tmp_L2;



endPoint4 = endPoint3 - length3 - u;

endx4 = startx3;
endy4 = starty3;


startx4 = endx4 + 2*length4*tan(var_divreq_x*DEG2RAD);
starty4 = endy4 + 2*length4*tan(var_divreq_y*DEG2RAD);


endPoint5 = endPoint4 - length4 - u;

endx5 = startx4;
endy5 = starty4;


var_divreq_x = atan(endx5*tan(var_divreq_x*DEG2RAD)/startx5)*RAD2DEG;
var_divreq_y = atan(endy5*tan(var_divreq_y*DEG2RAD)/starty5)*RAD2DEG;



if (startx5 > endx5) {
smallaxis_x5 = 0.5*startx5+ 0.5*smallaxis_x_factor5*(max_smallaxis_x5-startx5);
}
else {
smallaxis_x5 = 0.5*endx5+ 0.5*smallaxis_x_factor5*(max_smallaxis_x5-endx5);
}
if (starty5 > endy5) {
smallaxis_y5 = 0.5*starty5+ 0.5*smallaxis_y_factor5*(max_smallaxis_y5-starty5);
}
else {
smallaxis_y5 = 0.5*endy5+ 0.5*smallaxis_y_factor5*(max_smallaxis_y5-endy5);
}


// calculating focus points for elliptic guide 5 x direction
tmp_w1=startx5;
tmp_w2=endx5;
tmp_b=smallaxis_x5;
tmp_L=length5;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Linx5=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Loutx5= Linx5+tmp_L1-tmp_L2;

// calculating focus points for elliptic guide 5 y direction
tmp_w1=starty5;
tmp_w2=endy5;
tmp_b=smallaxis_y5;
tmp_L=length5;

tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));
tmp_L1=tmp_L/(1+tmp_k);
tmp_L2=tmp_L-tmp_L1;
tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));
Liny5=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);
Louty5= Liny5+tmp_L1-tmp_L2;






// Temporary fake start point, for this specific los section
startXdirec[7][1] = 1.0;
startXdirec[7][2] = 0.0;
startYdirec[7][1] = 1.0;
startYdirec[7][2] = 0.0;

startXposition[7][1] = 0.0;
startXposition[7][2] = 0.0;
startYposition[7][1] = 0.0;
startYposition[7][2] = 0.0;

startXdirec[6][1] = startXdirec[7][1];
startXdirec[6][2] = startXdirec[7][2];
startYdirec[6][1] = startYdirec[7][1];
startYdirec[6][2] = startYdirec[7][2];
startXposition[6][1] = startXposition[7][1] + startXdirec[7][1]*0*length5;
startXposition[6][2] = startXposition[7][2] + startXdirec[7][2]*0*length5;
startYposition[6][1] = startYposition[7][1] + startYdirec[7][1]*0*length5;
startYposition[6][2] = startYposition[7][2] + startYdirec[7][2]*0*length5;

position = 0*length5;
focus_e = length5 + Loutx5;
elength = focus_e + Linx5;
width  = smallaxis_x5*sqrt(1-(((position+Linx5)-elength/2)/(elength/2))*(((position+Linx5)-elength/2)/(elength/2)));
focus_e = length5 + Louty5;
elength = focus_e + Liny5;
height = smallaxis_y5*sqrt(1-(((position+Liny5)-elength/2)/(elength/2))*(((position+Liny5)-elength/2)/(elength/2)));
startxpoint[6][1][1] = startXposition[6][1]+startXdirec[6][2]*width;
startxpoint[6][2][1] = startXposition[6][2]-startXdirec[6][1]*width;
startxpoint[6][1][2] = startXposition[6][1]-startXdirec[6][2]*width;
startxpoint[6][2][2] = startXposition[6][2]+startXdirec[6][1]*width;
startypoint[6][1][1] = startYposition[6][1]+startYdirec[6][2]*height;
startypoint[6][2][1] = startYposition[6][2]-startYdirec[6][1]*height;
startypoint[6][1][2] = startYposition[6][1]-startYdirec[6][2]*height;
startypoint[6][2][2] = startYposition[6][2]+startYdirec[6][1]*height;



startXdirec[5][1] = startXdirec[7][1];
startXdirec[5][2] = startXdirec[7][2];
startYdirec[5][1] = startYdirec[7][1];
startYdirec[5][2] = startYdirec[7][2];
startXposition[5][1] = startXposition[7][1] + startXdirec[7][1]*length5;
startXposition[5][2] = startXposition[7][2] + startXdirec[7][2]*length5;
startYposition[5][1] = startYposition[7][1] + startYdirec[7][1]*length5;
startYposition[5][2] = startYposition[7][2] + startYdirec[7][2]*length5;

startxpoint[5][1][1] = startXposition[5][1]+startXdirec[5][2]*0.5*startx4;
startxpoint[5][2][1] = startXposition[5][2]-startXdirec[5][1]*0.5*startx4;
startxpoint[5][1][2] = startXposition[5][1]-startXdirec[5][2]*0.5*startx4;
startxpoint[5][2][2] = startXposition[5][2]+startXdirec[5][1]*0.5*startx4;
startypoint[5][1][1] = startYposition[5][1]+startYdirec[5][2]*0.5*starty4;
startypoint[5][2][1] = startYposition[5][2]-startYdirec[5][1]*0.5*starty4;
startypoint[5][1][2] = startYposition[5][1]-startYdirec[5][2]*0.5*starty4;
startypoint[5][2][2] = startYposition[5][2]+startYdirec[5][1]*0.5*starty4;




startXdirec[4][1] = startXdirec[5][1];
startXdirec[4][2] = startXdirec[5][2];
startYdirec[4][1] = startYdirec[5][1];
startYdirec[4][2] = startYdirec[5][2];
startXposition[4][1] = startXposition[5][1] + startXdirec[5][1]*length4;
startXposition[4][2] = startXposition[5][2] + startXdirec[5][2]*length4;
startYposition[4][1] = startYposition[5][1] + startYdirec[5][1]*length4;
startYposition[4][2] = startYposition[5][2] + startYdirec[5][2]*length4;

startxpoint[4][1][1] = startXposition[4][1]+startXdirec[4][2]*0.5*startx3;
startxpoint[4][2][1] = startXposition[4][2]-startXdirec[4][1]*0.5*startx3;
startxpoint[4][1][2] = startXposition[4][1]-startXdirec[4][2]*0.5*startx3;
startxpoint[4][2][2] = startXposition[4][2]+startXdirec[4][1]*0.5*startx3;
startypoint[4][1][1] = startYposition[4][1]+startYdirec[4][2]*0.5*starty3;
startypoint[4][2][1] = startYposition[4][2]-startYdirec[4][1]*0.5*starty3;
startypoint[4][1][2] = startYposition[4][1]-startYdirec[4][2]*0.5*starty3;
startypoint[4][2][2] = startYposition[4][2]+startYdirec[4][1]*0.5*starty3;




startXdirec[3][1] = startXdirec[4][1];
startXdirec[3][2] = startXdirec[4][2];
startYdirec[3][1] = startYdirec[4][1];
startYdirec[3][2] = startYdirec[4][2];
startXposition[3][1] = startXposition[4][1] + startXdirec[4][1]*length3;
startXposition[3][2] = startXposition[4][2] + startXdirec[4][2]*length3;
startYposition[3][1] = startYposition[4][1] + startYdirec[4][1]*length3;
startYposition[3][2] = startYposition[4][2] + startYdirec[4][2]*length3;

startxpoint[3][1][1] = startXposition[3][1]+startXdirec[3][2]*0.5*startx2;
startxpoint[3][2][1] = startXposition[3][2]-startXdirec[3][1]*0.5*startx2;
startxpoint[3][1][2] = startXposition[3][1]-startXdirec[3][2]*0.5*startx2;
startxpoint[3][2][2] = startXposition[3][2]+startXdirec[3][1]*0.5*startx2;
startypoint[3][1][1] = startYposition[3][1]+startYdirec[3][2]*0.5*starty2;
startypoint[3][2][1] = startYposition[3][2]-startYdirec[3][1]*0.5*starty2;
startypoint[3][1][2] = startYposition[3][1]-startYdirec[3][2]*0.5*starty2;
startypoint[3][2][2] = startYposition[3][2]+startYdirec[3][1]*0.5*starty2;

curveXcenter2[1] = startXposition[3][1] + startXdirec[3][2] * curve_radius2;
curveXcenter2[2] = startXposition[3][2] - startXdirec[3][1] * curve_radius2;
curve_small_radius2 = curve_radius2 - 0.5*startx2;



rot_factor_g2 = 1;
cosrot3 = cos(rot2*DEG2RAD);
sinrot3 = sin(rot2*DEG2RAD);
DeltaA2 = curve_radius2*(1-cos(rot2*DEG2RAD));
DeltaB2 = curve_radius2*sin(rot2*DEG2RAD);
startXdirec[2][1] = startXdirec[3][1]*cosrot3 + startXdirec[3][2]*sinrot3;
startXdirec[2][2] = - startXdirec[3][1]*sinrot3 + startXdirec[3][2]*cosrot3;
startXposition[2][1] = startXposition[3][1] + startXdirec[3][2]*DeltaA2 + startXdirec[3][1]*DeltaB2;
startXposition[2][2] = startXposition[3][2] - startXdirec[3][1]*DeltaA2 + startXdirec[3][2]*DeltaB2;
startYdirec[2][1] = startYdirec[3][1];
startYdirec[2][2] = startYdirec[3][2];
startYposition[2][1] = startYposition[3][1] + startYdirec[3][2]*DeltaA2 + startYdirec[3][1]*DeltaB2;
startYposition[2][2] = startYposition[3][2] + startYdirec[3][2]*length2*rot_factor_g2;

startxpoint[2][1][1] = startXposition[2][1]+startXdirec[2][2]*0.5*startx1;
startxpoint[2][2][1] = startXposition[2][2]-startXdirec[2][1]*0.5*startx1;
startxpoint[2][1][2] = startXposition[2][1]-startXdirec[2][2]*0.5*startx1;
startxpoint[2][2][2] = startXposition[2][2]+startXdirec[2][1]*0.5*startx1;
startypoint[2][1][1] = startYposition[2][1]+startYdirec[2][2]*0.5*starty1;
startypoint[2][2][1] = startYposition[2][2]-startYdirec[2][1]*0.5*starty1;
startypoint[2][1][2] = startYposition[2][1]-startYdirec[2][2]*0.5*starty1;
startypoint[2][2][2] = startYposition[2][2]+startYdirec[2][1]*0.5*starty1;




startXdirec[1][1] = startXdirec[2][1];
startXdirec[1][2] = startXdirec[2][2];
startYdirec[1][1] = startYdirec[2][1];
startYdirec[1][2] = startYdirec[2][2];
startXposition[1][1] = startXposition[2][1] + startXdirec[2][1]*(length1 - 6);
startXposition[1][2] = startXposition[2][2] + startXdirec[2][2]*(length1 - 6);
startYposition[1][1] = startYposition[2][1] + startYdirec[2][1]*(length1 - 6);
startYposition[1][2] = startYposition[2][2] + startYdirec[2][2]*(length1 - 6);

position = length1 - 6;
focus_e = length1 + Loutx1;
elength = focus_e + Linx1;
width  = smallaxis_x1*sqrt(1-(((position+Linx1)-elength/2)/(elength/2))*(((position+Linx1)-elength/2)/(elength/2)));
focus_e = length1 + Louty1;
elength = focus_e + Liny1;
height = smallaxis_y1*sqrt(1-(((position+Liny1)-elength/2)/(elength/2))*(((position+Liny1)-elength/2)/(elength/2)));
startxpoint[1][1][1] = startXposition[1][1]+startXdirec[1][2]*width;
startxpoint[1][2][1] = startXposition[1][2]-startXdirec[1][1]*width;
startxpoint[1][1][2] = startXposition[1][1]-startXdirec[1][2]*width;
startxpoint[1][2][2] = startXposition[1][2]+startXdirec[1][1]*width;
startypoint[1][1][1] = startYposition[1][1]+startYdirec[1][2]*height;
startypoint[1][2][1] = startYposition[1][2]-startYdirec[1][1]*height;
startypoint[1][1][2] = startYposition[1][1]-startYdirec[1][2]*height;
startypoint[1][2][2] = startYposition[1][2]+startYdirec[1][1]*height;
startXdirec[0][1] = startXdirec[2][1];
startXdirec[0][2] = startXdirec[2][2];
startYdirec[0][1] = startYdirec[2][1];
startYdirec[0][2] = startYdirec[2][2];
startXposition[0][1] = startXposition[2][1] + startXdirec[2][1]*length1;
startXposition[0][2] = startXposition[2][2] + startXdirec[2][2]*length1;
startYposition[0][1] = startYposition[2][1] + startYdirec[2][1]*length1;
startYposition[0][2] = startYposition[2][2] + startYdirec[2][2]*length1;
startxpoint[0][1][1] = startXposition[0][1]+startXdirec[0][2]*0.5*endx1;
startxpoint[0][2][1] = startXposition[0][2]-startXdirec[0][1]*0.5*endx1;
startxpoint[0][1][2] = startXposition[0][1]-startXdirec[0][2]*0.5*endx1;
startxpoint[0][2][2] = startXposition[0][2]+startXdirec[0][1]*0.5*endx1;
startypoint[0][1][1] = startYposition[0][1]+startYdirec[0][2]*0.5*endy1;
startypoint[0][2][1] = startYposition[0][2]-startYdirec[0][1]*0.5*endy1;
startypoint[0][1][2] = startYposition[0][1]-startYdirec[0][2]*0.5*endy1;
startypoint[0][2][2] = startYposition[0][2]+startYdirec[0][1]*0.5*endy1;

los_logic = 1; // assume line of sight
for (n1=1;n1<6;++n1) {
    for (n2=n1+1;n2<7;++n2) {
        Z1[1]=startxpoint[n1][1][1];
        X1[1]=startxpoint[n1][2][1];
        Z2[1]=startxpoint[n2][1][1];
        X2[1]=startxpoint[n2][2][1];
        Z1[2]=startxpoint[n1][1][2];
        X1[2]=startxpoint[n1][2][2];
        Z2[2]=startxpoint[n2][1][2];
        X2[2]=startxpoint[n2][2][2];
        Z1[3]=startxpoint[n1][1][1];
        X1[3]=startxpoint[n1][2][1];
        Z2[3]=startxpoint[n2][1][2];
        X2[3]=startxpoint[n2][2][2];
        Z1[4]=startxpoint[n1][1][2];
        X1[4]=startxpoint[n1][2][2];
        Z2[4]=startxpoint[n2][1][1];
        X2[4]=startxpoint[n2][2][1];
        for (line=1;line<5;++line) { // 4 lines need to be checked
            tmp_double=(Z1[line]-Z2[line])/100;
            a[line]=(X1[line]-X2[line])/tmp_double;
            tmp_double=a[line]*Z1[line];
            b[line]=X1[line]-tmp_double/100;
            los_tmp[line]=1;
            for (n_check=1;n_check<7;++n_check) {
                // check if the line goes between the points
                if (n_check != n1 && n_check != n2 && ((a[line]*startxpoint[n_check][1][1])/100+b[line]<startxpoint[n_check][2][1] || (a[line]*startxpoint[n_check][1][2])/100+b[line]>startxpoint[n_check][2][2])) {
                    los_tmp[line]=0; // line of sight blocked for line l!
                }
            }
            // x1 and x2 chosen as start and end of curved section
            x1_circ = startXposition[3][1] - curveXcenter2[1];
            x2_circ = startXposition[2][1] - curveXcenter2[1];
            y1_circ = a[line]/100*startXposition[3][1] + b[line] - curveXcenter2[2];
            y2_circ = a[line]/100*startXposition[2][1] + b[line] - curveXcenter2[2];
            dx_circ = x2_circ - x1_circ;
            dy_circ = y2_circ - y1_circ;
            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);
            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;
            if (dy_circ >= 0)
              sign_dy = 1;
            else
              sign_dy = -1;
            k_circ = curve_small_radius2*curve_small_radius2* dr_circ * dr_circ - D_circ*D_circ;

            if (k_circ > 0){
             x_solution[0] = (D_circ*dy_circ + sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);
             x_solution[1] = (D_circ*dy_circ - sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);
             y_solution[0] = (-D_circ*dx_circ + fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);
             y_solution[1] = (-D_circ*dx_circ - fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);

             for (ii=0;ii<2;++ii) {
              if (x_solution[ii] > x1_circ && x_solution[ii] < x2_circ) {
                  los_tmp[line]=0; // line of sight blocked for line l!
              }
             }
            }
        }
        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {
            los_logic_single[n1][n2]=0;
        }
        else {
            los_logic_single[n1][n2]=1;
        }
    }
}
los_check=0;
for (n1=1;n1<6;++n1) {
    for (n2=n1+1;n2<7;++n2) {
        los_check=los_check+los_logic_single[n1][n2];
    }
}
if (los_check==0) {
    los_logic = 0;
}
}







startXdirec[7][1] = 1.0;
startXdirec[7][2] = 0.0;
startYdirec[7][1] = 1.0;
startYdirec[7][2] = 0.0;

startXposition[7][1] = endPoint5 - length5;
startXposition[7][2] = 0.0;
startYposition[7][1] = endPoint5 - length5;
startYposition[7][2] = 0.0;

startxpoint[7][1][1] = startXposition[7][1]+startXdirec[7][2]*0.5*startx5;
startxpoint[7][2][1] = startXposition[7][2]-startXdirec[7][1]*0.5*startx5;
startxpoint[7][1][2] = startXposition[7][1]-startXdirec[7][2]*0.5*startx5;
startxpoint[7][2][2] = startXposition[7][2]+startXdirec[7][1]*0.5*startx5;
startypoint[7][1][1] = startYposition[7][1]+startYdirec[7][2]*0.5*starty5;
startypoint[7][2][1] = startYposition[7][2]-startYdirec[7][1]*0.5*starty5;
startypoint[7][1][2] = startYposition[7][1]-startYdirec[7][2]*0.5*starty5;
startypoint[7][2][2] = startYposition[7][2]+startYdirec[7][1]*0.5*starty5;




startXdirec[6][1] = startXdirec[7][1];
startXdirec[6][2] = startXdirec[7][2];
startYdirec[6][1] = startYdirec[7][1];
startYdirec[6][2] = startYdirec[7][2];
startXposition[6][1] = startXposition[7][1] + startXdirec[7][1]*0*length5;
startXposition[6][2] = startXposition[7][2] + startXdirec[7][2]*0*length5;
startYposition[6][1] = startYposition[7][1] + startYdirec[7][1]*0*length5;
startYposition[6][2] = startYposition[7][2] + startYdirec[7][2]*0*length5;

position = 0*length5;
focus_e = length5 + Loutx5;
elength = focus_e + Linx5;
width  = smallaxis_x5*sqrt(1-(((position+Linx5)-elength/2)/(elength/2))*(((position+Linx5)-elength/2)/(elength/2)));
focus_e = length5 + Louty5;
elength = focus_e + Liny5;
height = smallaxis_y5*sqrt(1-(((position+Liny5)-elength/2)/(elength/2))*(((position+Liny5)-elength/2)/(elength/2)));
startxpoint[6][1][1] = startXposition[6][1]+startXdirec[6][2]*width;
startxpoint[6][2][1] = startXposition[6][2]-startXdirec[6][1]*width;
startxpoint[6][1][2] = startXposition[6][1]-startXdirec[6][2]*width;
startxpoint[6][2][2] = startXposition[6][2]+startXdirec[6][1]*width;
startypoint[6][1][1] = startYposition[6][1]+startYdirec[6][2]*height;
startypoint[6][2][1] = startYposition[6][2]-startYdirec[6][1]*height;
startypoint[6][1][2] = startYposition[6][1]-startYdirec[6][2]*height;
startypoint[6][2][2] = startYposition[6][2]+startYdirec[6][1]*height;



startXdirec[5][1] = startXdirec[7][1];
startXdirec[5][2] = startXdirec[7][2];
startYdirec[5][1] = startYdirec[7][1];
startYdirec[5][2] = startYdirec[7][2];
startXposition[5][1] = startXposition[7][1] + startXdirec[7][1]*length5;
startXposition[5][2] = startXposition[7][2] + startXdirec[7][2]*length5;
startYposition[5][1] = startYposition[7][1] + startYdirec[7][1]*length5;
startYposition[5][2] = startYposition[7][2] + startYdirec[7][2]*length5;

startxpoint[5][1][1] = startXposition[5][1]+startXdirec[5][2]*0.5*startx4;
startxpoint[5][2][1] = startXposition[5][2]-startXdirec[5][1]*0.5*startx4;
startxpoint[5][1][2] = startXposition[5][1]-startXdirec[5][2]*0.5*startx4;
startxpoint[5][2][2] = startXposition[5][2]+startXdirec[5][1]*0.5*startx4;
startypoint[5][1][1] = startYposition[5][1]+startYdirec[5][2]*0.5*starty4;
startypoint[5][2][1] = startYposition[5][2]-startYdirec[5][1]*0.5*starty4;
startypoint[5][1][2] = startYposition[5][1]-startYdirec[5][2]*0.5*starty4;
startypoint[5][2][2] = startYposition[5][2]+startYdirec[5][1]*0.5*starty4;




startXdirec[4][1] = startXdirec[5][1];
startXdirec[4][2] = startXdirec[5][2];
startYdirec[4][1] = startYdirec[5][1];
startYdirec[4][2] = startYdirec[5][2];
startXposition[4][1] = startXposition[5][1] + startXdirec[5][1]*length4;
startXposition[4][2] = startXposition[5][2] + startXdirec[5][2]*length4;
startYposition[4][1] = startYposition[5][1] + startYdirec[5][1]*length4;
startYposition[4][2] = startYposition[5][2] + startYdirec[5][2]*length4;

startxpoint[4][1][1] = startXposition[4][1]+startXdirec[4][2]*0.5*startx3;
startxpoint[4][2][1] = startXposition[4][2]-startXdirec[4][1]*0.5*startx3;
startxpoint[4][1][2] = startXposition[4][1]-startXdirec[4][2]*0.5*startx3;
startxpoint[4][2][2] = startXposition[4][2]+startXdirec[4][1]*0.5*startx3;
startypoint[4][1][1] = startYposition[4][1]+startYdirec[4][2]*0.5*starty3;
startypoint[4][2][1] = startYposition[4][2]-startYdirec[4][1]*0.5*starty3;
startypoint[4][1][2] = startYposition[4][1]-startYdirec[4][2]*0.5*starty3;
startypoint[4][2][2] = startYposition[4][2]+startYdirec[4][1]*0.5*starty3;




startXdirec[3][1] = startXdirec[4][1];
startXdirec[3][2] = startXdirec[4][2];
startYdirec[3][1] = startYdirec[4][1];
startYdirec[3][2] = startYdirec[4][2];
startXposition[3][1] = startXposition[4][1] + startXdirec[4][1]*length3;
startXposition[3][2] = startXposition[4][2] + startXdirec[4][2]*length3;
startYposition[3][1] = startYposition[4][1] + startYdirec[4][1]*length3;
startYposition[3][2] = startYposition[4][2] + startYdirec[4][2]*length3;

startxpoint[3][1][1] = startXposition[3][1]+startXdirec[3][2]*0.5*startx2;
startxpoint[3][2][1] = startXposition[3][2]-startXdirec[3][1]*0.5*startx2;
startxpoint[3][1][2] = startXposition[3][1]-startXdirec[3][2]*0.5*startx2;
startxpoint[3][2][2] = startXposition[3][2]+startXdirec[3][1]*0.5*startx2;
startypoint[3][1][1] = startYposition[3][1]+startYdirec[3][2]*0.5*starty2;
startypoint[3][2][1] = startYposition[3][2]-startYdirec[3][1]*0.5*starty2;
startypoint[3][1][2] = startYposition[3][1]-startYdirec[3][2]*0.5*starty2;
startypoint[3][2][2] = startYposition[3][2]+startYdirec[3][1]*0.5*starty2;

curveXcenter2[1] = startXposition[3][1] + startXdirec[3][2] * curve_radius2;
curveXcenter2[2] = startXposition[3][2] - startXdirec[3][1] * curve_radius2;
curve_small_radius2 = curve_radius2 - 0.5*startx2;



rot_factor_g2 = 1;
cosrot3 = cos(rot2*DEG2RAD);
sinrot3 = sin(rot2*DEG2RAD);
DeltaA2 = curve_radius2*(1-cos(rot2*DEG2RAD));
DeltaB2 = curve_radius2*sin(rot2*DEG2RAD);
startXdirec[2][1] = startXdirec[3][1]*cosrot3 + startXdirec[3][2]*sinrot3;
startXdirec[2][2] = - startXdirec[3][1]*sinrot3 + startXdirec[3][2]*cosrot3;
startXposition[2][1] = startXposition[3][1] + startXdirec[3][2]*DeltaA2 + startXdirec[3][1]*DeltaB2;
startXposition[2][2] = startXposition[3][2] - startXdirec[3][1]*DeltaA2 + startXdirec[3][2]*DeltaB2;
startYdirec[2][1] = startYdirec[3][1];
startYdirec[2][2] = startYdirec[3][2];
startYposition[2][1] = startYposition[3][1] + startYdirec[3][2]*DeltaA2 + startYdirec[3][1]*DeltaB2;
startYposition[2][2] = startYposition[3][2] + startYdirec[3][2]*length2*rot_factor_g2;

startxpoint[2][1][1] = startXposition[2][1]+startXdirec[2][2]*0.5*startx1;
startxpoint[2][2][1] = startXposition[2][2]-startXdirec[2][1]*0.5*startx1;
startxpoint[2][1][2] = startXposition[2][1]-startXdirec[2][2]*0.5*startx1;
startxpoint[2][2][2] = startXposition[2][2]+startXdirec[2][1]*0.5*startx1;
startypoint[2][1][1] = startYposition[2][1]+startYdirec[2][2]*0.5*starty1;
startypoint[2][2][1] = startYposition[2][2]-startYdirec[2][1]*0.5*starty1;
startypoint[2][1][2] = startYposition[2][1]-startYdirec[2][2]*0.5*starty1;
startypoint[2][2][2] = startYposition[2][2]+startYdirec[2][1]*0.5*starty1;




startXdirec[1][1] = startXdirec[2][1];
startXdirec[1][2] = startXdirec[2][2];
startYdirec[1][1] = startYdirec[2][1];
startYdirec[1][2] = startYdirec[2][2];
startXposition[1][1] = startXposition[2][1] + startXdirec[2][1]*(length1 - 6);
startXposition[1][2] = startXposition[2][2] + startXdirec[2][2]*(length1 - 6);
startYposition[1][1] = startYposition[2][1] + startYdirec[2][1]*(length1 - 6);
startYposition[1][2] = startYposition[2][2] + startYdirec[2][2]*(length1 - 6);

position = length1 - 6;
focus_e = length1 + Loutx1;
elength = focus_e + Linx1;
width  = smallaxis_x1*sqrt(1-(((position+Linx1)-elength/2)/(elength/2))*(((position+Linx1)-elength/2)/(elength/2)));
focus_e = length1 + Louty1;
elength = focus_e + Liny1;
height = smallaxis_y1*sqrt(1-(((position+Liny1)-elength/2)/(elength/2))*(((position+Liny1)-elength/2)/(elength/2)));
startxpoint[1][1][1] = startXposition[1][1]+startXdirec[1][2]*width;
startxpoint[1][2][1] = startXposition[1][2]-startXdirec[1][1]*width;
startxpoint[1][1][2] = startXposition[1][1]-startXdirec[1][2]*width;
startxpoint[1][2][2] = startXposition[1][2]+startXdirec[1][1]*width;
startypoint[1][1][1] = startYposition[1][1]+startYdirec[1][2]*height;
startypoint[1][2][1] = startYposition[1][2]-startYdirec[1][1]*height;
startypoint[1][1][2] = startYposition[1][1]-startYdirec[1][2]*height;
startypoint[1][2][2] = startYposition[1][2]+startYdirec[1][1]*height;
startXdirec[0][1] = startXdirec[2][1];
startXdirec[0][2] = startXdirec[2][2];
startYdirec[0][1] = startYdirec[2][1];
startYdirec[0][2] = startYdirec[2][2];
startXposition[0][1] = startXposition[2][1] + startXdirec[2][1]*length1;
startXposition[0][2] = startXposition[2][2] + startXdirec[2][2]*length1;
startYposition[0][1] = startYposition[2][1] + startYdirec[2][1]*length1;
startYposition[0][2] = startYposition[2][2] + startYdirec[2][2]*length1;
startxpoint[0][1][1] = startXposition[0][1]+startXdirec[0][2]*0.5*endx1;
startxpoint[0][2][1] = startXposition[0][2]-startXdirec[0][1]*0.5*endx1;
startxpoint[0][1][2] = startXposition[0][1]-startXdirec[0][2]*0.5*endx1;
startxpoint[0][2][2] = startXposition[0][2]+startXdirec[0][1]*0.5*endx1;
startypoint[0][1][1] = startYposition[0][1]+startYdirec[0][2]*0.5*endy1;
startypoint[0][2][1] = startYposition[0][2]-startYdirec[0][1]*0.5*endy1;
startypoint[0][1][2] = startYposition[0][1]-startYdirec[0][2]*0.5*endy1;
startypoint[0][2][2] = startYposition[0][2]+startYdirec[0][1]*0.5*endy1;


fp = fopen(file_name,"w");

fprintf(fp,"E\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[7][1][1],startxpoint[7][2][1],startxpoint[7][1][2],startxpoint[7][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[6][1][1],startxpoint[6][2][1],startxpoint[6][1][2],startxpoint[6][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[7][1][1],startypoint[7][2][1],startypoint[7][1][2],startypoint[7][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[6][1][1],startypoint[6][2][1],startypoint[6][1][2],startypoint[6][2][2]);
fprintf(fp,"xpos: %lf \t %lf\n",startXposition[7][1],startXposition[7][2]);
fprintf(fp,"ypos: %lf \t %lf\n",startYposition[7][1],startYposition[7][2]);
fprintf(fp,"xdir: %lf \t %lf\n",startXdirec[7][1],startXdirec[7][2]);
fprintf(fp,"ydir: %lf \t %lf\n",startYdirec[7][1],startYdirec[7][2]);
fprintf(fp,"xpars: %lf \t %lf \t %lf \t %lf\n",length5,smallaxis_x5,Linx5,Loutx5);
fprintf(fp,"ypars: %lf \t %lf \t %lf \t %lf\n",length5,smallaxis_y5,Liny5,Louty5);

fprintf(fp,"E\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[6][1][1],startxpoint[6][2][1],startxpoint[6][1][2],startxpoint[6][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[5][1][1],startxpoint[5][2][1],startxpoint[5][1][2],startxpoint[5][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[6][1][1],startypoint[6][2][1],startypoint[6][1][2],startypoint[6][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[5][1][1],startypoint[5][2][1],startypoint[5][1][2],startypoint[5][2][2]);
fprintf(fp,"xpos: %lf \t %lf\n",startXposition[6][1],startXposition[6][2]);
fprintf(fp,"ypos: %lf \t %lf\n",startYposition[6][1],startYposition[6][2]);
fprintf(fp,"xdir: %lf \t %lf\n",startXdirec[6][1],startXdirec[6][2]);
fprintf(fp,"ydir: %lf \t %lf\n",startYdirec[6][1],startYdirec[6][2]);
fprintf(fp,"xpars: %lf \t %lf \t %lf \t %lf\n",length5,smallaxis_x5,Linx5,Loutx5);
fprintf(fp,"ypars: %lf \t %lf \t %lf \t %lf\n",length5,smallaxis_y5,Liny5,Louty5);
fprintf(fp,"Los s\n");

fprintf(fp,"G\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[5][1][1],startxpoint[5][2][1],startxpoint[5][1][2],startxpoint[5][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[4][1][1],startxpoint[4][2][1],startxpoint[4][1][2],startxpoint[4][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[5][1][1],startypoint[5][2][1],startypoint[5][1][2],startypoint[5][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[4][1][1],startypoint[4][2][1],startypoint[4][1][2],startypoint[4][2][2]);

fprintf(fp,"E\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[4][1][1],startxpoint[4][2][1],startxpoint[4][1][2],startxpoint[4][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[3][1][1],startxpoint[3][2][1],startxpoint[3][1][2],startxpoint[3][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[4][1][1],startypoint[4][2][1],startypoint[4][1][2],startypoint[4][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[3][1][1],startypoint[3][2][1],startypoint[3][1][2],startypoint[3][2][2]);
fprintf(fp,"xpos: %lf \t %lf\n",startXposition[4][1],startXposition[4][2]);
fprintf(fp,"ypos: %lf \t %lf\n",startYposition[4][1],startYposition[4][2]);
fprintf(fp,"xdir: %lf \t %lf\n",startXdirec[4][1],startXdirec[4][2]);
fprintf(fp,"ydir: %lf \t %lf\n",startYdirec[4][1],startYdirec[4][2]);
fprintf(fp,"xpars: %lf \t %lf \t %lf \t %lf\n",length3,smallaxis_x3,Linx3,Loutx3);
fprintf(fp,"ypars: %lf \t %lf \t %lf \t %lf\n",length3,smallaxis_y3,Liny3,Louty3);

fprintf(fp,"C\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[3][1][1],startxpoint[3][2][1],startxpoint[3][1][2],startxpoint[3][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[2][1][1],startxpoint[2][2][1],startxpoint[2][1][2],startxpoint[2][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[3][1][1],startypoint[3][2][1],startypoint[3][1][2],startypoint[3][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[2][1][1],startypoint[2][2][1],startypoint[2][1][2],startypoint[2][2][2]);
fprintf(fp,"x: %lf \t %lf \t %lf \t %lf \t %lf\n",curve_radius2,curve_small_radius2,curveXcenter2[1],curveXcenter2[2],channels2);

fprintf(fp,"E\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[2][1][1],startxpoint[2][2][1],startxpoint[2][1][2],startxpoint[2][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[1][1][1],startxpoint[1][2][1],startxpoint[1][1][2],startxpoint[1][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[2][1][1],startypoint[2][2][1],startypoint[2][1][2],startypoint[2][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[1][1][1],startypoint[1][2][1],startypoint[1][1][2],startypoint[1][2][2]);
fprintf(fp,"xpos: %lf \t %lf\n",startXposition[2][1],startXposition[2][2]);
fprintf(fp,"ypos: %lf \t %lf\n",startYposition[2][1],startYposition[2][2]);
fprintf(fp,"xdir: %lf \t %lf\n",startXdirec[2][1],startXdirec[2][2]);
fprintf(fp,"ydir: %lf \t %lf\n",startYdirec[2][1],startYdirec[2][2]);
fprintf(fp,"xpars: %lf \t %lf \t %lf \t %lf\n",length1,smallaxis_x1,Linx1,Loutx1);
fprintf(fp,"ypars: %lf \t %lf \t %lf \t %lf\n",length1,smallaxis_y1,Liny1,Louty1);

fprintf(fp,"E\n");
fprintf(fp,"x1: %lf \t %lf \t %lf \t %lf\n",startxpoint[1][1][1],startxpoint[1][2][1],startxpoint[1][1][2],startxpoint[1][2][2]);
fprintf(fp,"x2: %lf \t %lf \t %lf \t %lf\n",startxpoint[0][1][1],startxpoint[0][2][1],startxpoint[0][1][2],startxpoint[0][2][2]);
fprintf(fp,"y1: %lf \t %lf \t %lf \t %lf\n",startypoint[1][1][1],startypoint[1][2][1],startypoint[1][1][2],startypoint[1][2][2]);
fprintf(fp,"y2: %lf \t %lf \t %lf \t %lf\n",startypoint[0][1][1],startypoint[0][2][1],startypoint[0][1][2],startypoint[0][2][2]);
fprintf(fp,"xpos: %lf \t %lf\n",startXposition[1][1],startXposition[1][2]);
fprintf(fp,"ypos: %lf \t %lf\n",startYposition[1][1],startYposition[1][2]);
fprintf(fp,"xdir: %lf \t %lf\n",startXdirec[1][1],startXdirec[1][2]);
fprintf(fp,"ydir: %lf \t %lf\n",startYdirec[1][1],startYdirec[1][2]);
fprintf(fp,"xpars: %lf \t %lf \t %lf \t %lf\n",length1,smallaxis_x1,Linx1,Loutx1);
fprintf(fp,"ypars: %lf \t %lf \t %lf \t %lf\n",length1,smallaxis_y1,Liny1,Louty1);
fprintf(fp,"Los e\n");
fprintf(fp,"Sample\n");
fprintf(fp,"%lf \t %lf \t %lf\n",sizeX,sizeY,sample_dist);
fprintf(fp,"xp: %lf \t %lf\n",startXposition[0][1],startXposition[0][2]);
fprintf(fp,"yp: %lf \t %lf\n",startYposition[0][1],startYposition[0][2]);
fprintf(fp,"xd: %lf \t %lf\n",startXdirec[0][1],startXdirec[0][2]);
fprintf(fp,"yd: %lf \t %lf\n",startYdirec[0][1],startYdirec[0][2]);
fprintf(fp,"Moderator\n");
fprintf(fp,"%lf \t %lf\n",mod_x,mod_y);
fclose(fp);

double BE_THICK = 0.12; // Originally from Secondary Spectrometer
int COARSE = 1;
// Filter length*/
OUTFILTER_THICK = BE_THICK/0.05;


// Coarse collimator 
if (COARSE) {coarse = 1;}
else {coarse = 0;}


//Pixel start
pixelstart[0]=96;
pixelstart[1]=205;
pixelstart[2]=312;
pixelstart[3]=418;
pixelstart[4]=525;
pixelstart[5]=632;
pixelstart[6]=741;
pixelstart[7]=851;

pixelstop[0]=178;
pixelstop[1]=285;
pixelstop[2]=392;
pixelstop[3]=498;
pixelstop[4]=605;
pixelstop[5]=713;
pixelstop[6]=823;
pixelstop[7]=935;


// Analyzer positions
z_an1[1]=0.9115;
z_an1[2]=0.9209;
z_an1[3]=0.9300;
z_an1[4]=0.9388;
z_an1[5]=0.9474;
y_an1[1]=-0.0275;
y_an1[2]=-0.0140;
y_an1[3]=0.0000;
y_an1[4]=0.0145;
y_an1[5]=0.0297;
thet_an1[1]=47.17;
thet_an1[2]=48.03;
thet_an1[3]=48.90;
thet_an1[4]=49.79;
thet_an1[5]=50.69;
ana_l[1]=0.0720;
Ana1=0;

z_an2[1]=0.9744;
z_an2[2]=0.9843;
z_an2[3]=0.9939;
z_an2[4]=1.0033;
z_an2[5]=1.0123;
y_an2[1]=-0.0290;
y_an2[2]=-0.0148;
y_an2[3]=0.0000;
y_an2[4]=0.0153;
y_an2[5]=0.0312;
thet_an2[1]=45.51;
thet_an2[2]=46.35;
thet_an2[3]=47.21;
thet_an2[4]=48.09;
thet_an2[5]=48.98;
ana_l[2]=0.082;
Ana2=0;

z_an3[1]=1.0354;
z_an3[2]=1.0462;
z_an3[3]=1.0569;
z_an3[4]=1.0671;
z_an3[5]=1.0771;
y_an3[1]=-0.0311;
y_an3[2]=-0.0160;
y_an3[3]=0.0000;
y_an3[4]=0.0165;
y_an3[5]=0.0338;
thet_an3[1]=43.81;
thet_an3[2]=44.65;
thet_an3[3]=45.53;
thet_an3[4]=46.41;
thet_an3[5]=47.32;
ana_l[3]=0.092;
Ana3=0;

z_an4[1]=1.0959;
z_an4[2]=1.1078;
z_an4[3]=1.1195;
z_an4[4]=1.1309;
z_an4[5]=1.1420;
y_an4[1]=-0.0332;
y_an4[2]=-0.0170;
y_an4[3]=0.0000;
y_an4[4]=0.0177;
y_an4[5]=0.0362;
thet_an4[1]=42.11;
thet_an4[2]=42.96;
thet_an4[3]=43.84;
thet_an4[4]=44.74;
thet_an4[5]=45.66;
ana_l[4]=0.103;
Ana4=0;

z_an5[1]=1.1567;
z_an5[2]=1.1697;
z_an5[3]=1.1827;
z_an5[4]=1.1953;
z_an5[5]=1.2076;
y_an5[1]=-0.0350;
y_an5[2]=-0.0180;
y_an5[3]=0.0000;
y_an5[4]=0.0187;
y_an5[5]=0.0382;
thet_an5[1]=40.42;
thet_an5[2]=41.27;
thet_an5[3]=42.16;
thet_an5[4]=43.05;
thet_an5[5]=43.97;
ana_l[5]=0.113;
Ana5=0;

z_an6[1]=1.2170;
z_an6[2]=1.2313;
z_an6[3]=1.2456;
z_an6[4]=1.2595;
z_an6[5]=1.2732;
y_an6[1]=-0.0365;
y_an6[2]=-0.0188;
y_an6[3]=0.0000;
y_an6[4]=0.0195;
y_an6[5]=0.0402;
thet_an6[1]=38.75;
thet_an6[2]=39.60;
thet_an6[3]=40.47;
thet_an6[4]=41.36;
thet_an6[5]=42.28;
ana_l[6]=0.120;
Ana6=0;

z_an7[1]=1.2775;
z_an7[2]=1.2936;
z_an7[3]=1.3098;
z_an7[4]=1.3256;
z_an7[5]=1.3413;
y_an7[1]=-0.0389;
y_an7[2]=-0.0201;
y_an7[3]=0.0000;
y_an7[4]=0.0209;
y_an7[5]=0.0431;
thet_an7[1]=37.04;
thet_an7[2]=37.90;
thet_an7[3]=38.79;
thet_an7[4]=39.69;
thet_an7[5]=40.63;
ana_l[7]=0.129;
Ana7=0;

z_an8[1]=1.3385;
z_an8[2]=1.3565;
z_an8[3]=1.3747;
z_an8[4]=1.3925;
z_an8[5]=1.4102;
y_an8[1]=-0.0409;
y_an8[2]=-0.0212;
y_an8[3]=0.0000;
y_an8[4]=0.0221;
y_an8[5]=0.0456;
thet_an8[1]=35.35;
thet_an8[2]=36.21;
thet_an8[3]=37.10;
thet_an8[4]=38.01;
thet_an8[5]=38.95;
ana_l[8]=0.140;
Ana8=0;

ana_l[1]=ana_l[1]+ana_l_add;
ana_l[2]=ana_l[2]+ana_l_add;
ana_l[3]=ana_l[3]+ana_l_add;
ana_l[4]=ana_l[4]+ana_l_add;
ana_l[5]=ana_l[5]+ana_l_add;
ana_l[6]=ana_l[6]+ana_l_add;
ana_l[7]=ana_l[7]+ana_l_add;
ana_l[8]=ana_l[8]+ana_l_add;


// Detector positions rotation given in degrees

ang_1[0]=-3.33;
ang_1[1]=-2.22;
ang_1[2]=-1.11;
ang_1[3]=0;
ang_1[4]=1.11;
ang_1[5]=2.22;
ang_1[6]=3.33;

ang_2[0]=-2.775;
ang_2[1]=-1.665;
ang_2[2]=-0.555;
ang_2[3]=0.555;
ang_2[4]=1.665;
ang_2[5]=2.775;

dety_1=0.7; //
dety_2=0.71; //
det_cen=DET_cen;

A3_val = A3;//*3.141592/180.0;

A4_val = A4;//*3.141592/180.0;
//sourceDivergence = source_DivH;

for (i=0; i<7; i++){
	detx_1[i] = -sin(ang_1[i]*PI/180)*det_cen;
	detz_1[i] = cos(ang_1[i]*PI/180)*det_cen;
}

for (i=0; i<6; i++){
	detx_2[i] = -sin(ang_2[i]*PI/180)*det_cen;
	detz_2[i] = cos(ang_2[i]*PI/180)*det_cen;
}


RV_mono = 2*1.6*sin(DEG2RAD*thetaA1)*RV_mono_Bool;
RH_mono = 2*1.6/sin(DEG2RAD*thetaA1)*RH_mono_Bool;

}
#line 10517 "./FullInstrument_v4.7_frontend.c"
#undef samplePosY
#undef sampleRadius
#undef sampleHeight
#undef old_mono
#undef Be_filter
#undef RH_mono_Bool
#undef RV_mono_Bool
#undef file_name
#undef SAMPLE
#undef A4
#undef A3
#undef EI
#undef DeltaSourceE
#undef SourceE
#undef virtual_slit_d
#undef slitwidth
#undef mcposaPSI_CAMEA
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component Origin. */
  /* Setting parameters for component Origin. */
  SIG_MESSAGE("Origin (Init:SetPar)");
#line 39 "FullInstrument_v4.7_frontend.instr"
  if("NULL") strncpy(mccOrigin_profile, "NULL" ? "NULL" : "", 16384); else mccOrigin_profile[0]='\0';
#line 39 "FullInstrument_v4.7_frontend.instr"
  mccOrigin_percent = 10;
#line 39 "FullInstrument_v4.7_frontend.instr"
  mccOrigin_flag_save = 0;
#line 39 "FullInstrument_v4.7_frontend.instr"
  mccOrigin_minutes = 0;
#line 10559 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10566 "./FullInstrument_v4.7_frontend.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 1332 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1332 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1332 "FullInstrument_v4.7_frontend.instr"
    0);
#line 10575 "./FullInstrument_v4.7_frontend.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");
#line 120 "FullInstrument_v4.7_frontend.instr"
  mccsource_radius = 0.0;
#line 1339 "FullInstrument_v4.7_frontend.instr"
  mccsource_dist = 1.465;
#line 1337 "FullInstrument_v4.7_frontend.instr"
  mccsource_xw = 0.035;
#line 1337 "FullInstrument_v4.7_frontend.instr"
  mccsource_yh = 0.12;
#line 120 "FullInstrument_v4.7_frontend.instr"
  mccsource_E0 = 0;
#line 120 "FullInstrument_v4.7_frontend.instr"
  mccsource_dE = 0;
#line 120 "FullInstrument_v4.7_frontend.instr"
  mccsource_Lambda0 = 0;
#line 120 "FullInstrument_v4.7_frontend.instr"
  mccsource_dLambda = 0;
#line 1341 "FullInstrument_v4.7_frontend.instr"
  mccsource_I1 = ( 1.27e13 / 4 / PI );
#line 1337 "FullInstrument_v4.7_frontend.instr"
  mccsource_h = 0.135;
#line 1337 "FullInstrument_v4.7_frontend.instr"
  mccsource_w = 0.08;
#line 121 "FullInstrument_v4.7_frontend.instr"
  mccsource_gaussian = 0;
#line 121 "FullInstrument_v4.7_frontend.instr"
  mccsource_verbose = 0;
#line 1341 "FullInstrument_v4.7_frontend.instr"
  mccsource_T1 = 301.287;
#line 122 "FullInstrument_v4.7_frontend.instr"
  mccsource_flux_file_perAA = 0;
#line 122 "FullInstrument_v4.7_frontend.instr"
  mccsource_flux_file_log = 0;
#line 1339 "FullInstrument_v4.7_frontend.instr"
  mccsource_Lmin = lmin;
#line 1339 "FullInstrument_v4.7_frontend.instr"
  mccsource_Lmax = lmax;
#line 124 "FullInstrument_v4.7_frontend.instr"
  mccsource_Emin = 0;
#line 124 "FullInstrument_v4.7_frontend.instr"
  mccsource_Emax = 0;
#line 1342 "FullInstrument_v4.7_frontend.instr"
  mccsource_T2 = 105.655;
#line 1342 "FullInstrument_v4.7_frontend.instr"
  mccsource_I2 = ( 3.818e12 / 4 / PI );
#line 1343 "FullInstrument_v4.7_frontend.instr"
  mccsource_T3 = 25.379;
#line 1343 "FullInstrument_v4.7_frontend.instr"
  mccsource_I3 = ( 2.331e12 / 4 / PI );
#line 124 "FullInstrument_v4.7_frontend.instr"
  mccsource_length = 0;
#line 124 "FullInstrument_v4.7_frontend.instr"
  mccsource_phi_init = 0;
#line 124 "FullInstrument_v4.7_frontend.instr"
  mccsource_theta_init = 0;
#line 1344 "FullInstrument_v4.7_frontend.instr"
  mccsource_HEtailA = 8.306e11 / 4 / PI;
#line 1344 "FullInstrument_v4.7_frontend.instr"
  mccsource_HEtailL0 = -0.398;
#line 10644 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10651 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotasource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 1345 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1345 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1345 "FullInstrument_v4.7_frontend.instr"
    0);
#line 10662 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[2] = mcposasource;
  mccomp_posr[2] = mcposrsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component StartOfGuide. */
  /* Setting parameters for component StartOfGuide. */
  SIG_MESSAGE("StartOfGuide (Init:SetPar)");

  SIG_MESSAGE("StartOfGuide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10682 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotasource, mcrotaStartOfGuide);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotaStartOfGuide, mctr1, mcrotrStartOfGuide);
  mctc1 = coords_set(
#line 1349 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1349 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1349 "FullInstrument_v4.7_frontend.instr"
    guide_start);
#line 10693 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaStartOfGuide = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposasource, mcposaStartOfGuide);
  mcposrStartOfGuide = rot_apply(mcrotaStartOfGuide, mctc1);
  mcDEBUG_COMPONENT("StartOfGuide", mcposaStartOfGuide, mcrotaStartOfGuide)
  mccomp_posa[3] = mcposaStartOfGuide;
  mccomp_posr[3] = mcposrStartOfGuide;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component elliptical_guide_gravity5. */
  /* Setting parameters for component elliptical_guide_gravity5. */
  SIG_MESSAGE("elliptical_guide_gravity5 (Init:SetPar)");
#line 1355 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_l = length5;
#line 1356 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_xwidth = 2 * smallaxis_x5;
#line 1356 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_yheight = 2 * smallaxis_y5;
#line 1355 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_linxw = Linx5;
#line 1355 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_loutxw = Loutx5;
#line 1355 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_linyh = Liny5;
#line 1355 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_loutyh = Louty5;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_majorAxisxw = 0;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_minorAxisxw = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_majorAxisyh = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_minorAxisyh = 0;
#line 167 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_majorAxisoffsetxw = 0;
#line 168 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_majorAxisoffsetyh = 0;
#line 1356 "FullInstrument_v4.7_frontend.instr"
  if("mid") strncpy(mccelliptical_guide_gravity5_dimensionsAt, "mid" ? "mid" : "", 16384); else mccelliptical_guide_gravity5_dimensionsAt[0]='\0';
#line 170 "FullInstrument_v4.7_frontend.instr"
  if("ellipse") strncpy(mccelliptical_guide_gravity5_option, "ellipse" ? "ellipse" : "", 16384); else mccelliptical_guide_gravity5_option[0]='\0';
#line 1357 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_R0 = R05;
#line 1357 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_Qc = Qc5;
#line 1357 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_alpha = m5 * 0.24 + 2.1;
#line 1357 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_m = m5;
#line 1357 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_W = W5;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_alpharight = -1;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_mright = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_alphaleft = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_mleft = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_alphatop = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_mtop = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_alphabottom = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_mbottom = -1;
#line 176 "FullInstrument_v4.7_frontend.instr"
  if("on") strncpy(mccelliptical_guide_gravity5_verbose, "on" ? "on" : "", 16384); else mccelliptical_guide_gravity5_verbose[0]='\0';
#line 177 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_enableGravity = 1.0;
#line 178 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity5_curvature = 0;
#line 10769 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("elliptical_guide_gravity5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1359 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1359 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1359 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 10779 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaStartOfGuide, mcrotaelliptical_guide_gravity5);
  rot_transpose(mcrotaStartOfGuide, mctr1);
  rot_mul(mcrotaelliptical_guide_gravity5, mctr1, mcrotrelliptical_guide_gravity5);
  mctc1 = coords_set(
#line 1358 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1358 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1358 "FullInstrument_v4.7_frontend.instr"
    u);
#line 10790 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaStartOfGuide, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaelliptical_guide_gravity5 = coords_add(mcposaStartOfGuide, mctc2);
  mctc1 = coords_sub(mcposaStartOfGuide, mcposaelliptical_guide_gravity5);
  mcposrelliptical_guide_gravity5 = rot_apply(mcrotaelliptical_guide_gravity5, mctc1);
  mcDEBUG_COMPONENT("elliptical_guide_gravity5", mcposaelliptical_guide_gravity5, mcrotaelliptical_guide_gravity5)
  mccomp_posa[4] = mcposaelliptical_guide_gravity5;
  mccomp_posr[4] = mcposrelliptical_guide_gravity5;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component EndOfelement_5. */
  /* Setting parameters for component EndOfelement_5. */
  SIG_MESSAGE("EndOfelement_5 (Init:SetPar)");

  SIG_MESSAGE("EndOfelement_5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10810 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaelliptical_guide_gravity5, mcrotaEndOfelement_5);
  rot_transpose(mcrotaelliptical_guide_gravity5, mctr1);
  rot_mul(mcrotaEndOfelement_5, mctr1, mcrotrEndOfelement_5);
  mctc1 = coords_set(
#line 1362 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1362 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1362 "FullInstrument_v4.7_frontend.instr"
    length5 + 2 * u);
#line 10821 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaelliptical_guide_gravity5, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEndOfelement_5 = coords_add(mcposaelliptical_guide_gravity5, mctc2);
  mctc1 = coords_sub(mcposaelliptical_guide_gravity5, mcposaEndOfelement_5);
  mcposrEndOfelement_5 = rot_apply(mcrotaEndOfelement_5, mctc1);
  mcDEBUG_COMPONENT("EndOfelement_5", mcposaEndOfelement_5, mcrotaEndOfelement_5)
  mccomp_posa[5] = mcposaEndOfelement_5;
  mccomp_posr[5] = mcposrEndOfelement_5;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component EndOfelement_4. */
  /* Setting parameters for component EndOfelement_4. */
  SIG_MESSAGE("EndOfelement_4 (Init:SetPar)");

  SIG_MESSAGE("EndOfelement_4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10841 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_5, mcrotaEndOfelement_4);
  rot_transpose(mcrotaEndOfelement_5, mctr1);
  rot_mul(mcrotaEndOfelement_4, mctr1, mcrotrEndOfelement_4);
  mctc1 = coords_set(
#line 1367 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1367 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1367 "FullInstrument_v4.7_frontend.instr"
    length4);
#line 10852 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_5, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEndOfelement_4 = coords_add(mcposaEndOfelement_5, mctc2);
  mctc1 = coords_sub(mcposaEndOfelement_5, mcposaEndOfelement_4);
  mcposrEndOfelement_4 = rot_apply(mcrotaEndOfelement_4, mctc1);
  mcDEBUG_COMPONENT("EndOfelement_4", mcposaEndOfelement_4, mcrotaEndOfelement_4)
  mccomp_posa[6] = mcposaEndOfelement_4;
  mccomp_posr[6] = mcposrEndOfelement_4;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component elliptical_guide_gravity3. */
  /* Setting parameters for component elliptical_guide_gravity3. */
  SIG_MESSAGE("elliptical_guide_gravity3 (Init:SetPar)");
#line 1371 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_l = length3;
#line 1372 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_xwidth = 2 * smallaxis_x3;
#line 1372 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_yheight = 2 * smallaxis_y3;
#line 1371 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_linxw = Linx3;
#line 1371 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_loutxw = Loutx3;
#line 1371 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_linyh = Liny3;
#line 1371 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_loutyh = Louty3;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_majorAxisxw = 0;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_minorAxisxw = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_majorAxisyh = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_minorAxisyh = 0;
#line 167 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_majorAxisoffsetxw = 0;
#line 168 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_majorAxisoffsetyh = 0;
#line 1372 "FullInstrument_v4.7_frontend.instr"
  if("mid") strncpy(mccelliptical_guide_gravity3_dimensionsAt, "mid" ? "mid" : "", 16384); else mccelliptical_guide_gravity3_dimensionsAt[0]='\0';
#line 170 "FullInstrument_v4.7_frontend.instr"
  if("ellipse") strncpy(mccelliptical_guide_gravity3_option, "ellipse" ? "ellipse" : "", 16384); else mccelliptical_guide_gravity3_option[0]='\0';
#line 1373 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_R0 = R03;
#line 1373 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_Qc = Qc3;
#line 1373 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_alpha = m3 * 0.24 + 2.1;
#line 1373 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_m = m3;
#line 1373 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_W = W3;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_alpharight = -1;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_mright = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_alphaleft = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_mleft = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_alphatop = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_mtop = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_alphabottom = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_mbottom = -1;
#line 176 "FullInstrument_v4.7_frontend.instr"
  if("on") strncpy(mccelliptical_guide_gravity3_verbose, "on" ? "on" : "", 16384); else mccelliptical_guide_gravity3_verbose[0]='\0';
#line 177 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_enableGravity = 1.0;
#line 178 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity3_curvature = 0;
#line 10928 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("elliptical_guide_gravity3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1375 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1375 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1375 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 10938 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_4, mcrotaelliptical_guide_gravity3);
  rot_transpose(mcrotaEndOfelement_4, mctr1);
  rot_mul(mcrotaelliptical_guide_gravity3, mctr1, mcrotrelliptical_guide_gravity3);
  mctc1 = coords_set(
#line 1374 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1374 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1374 "FullInstrument_v4.7_frontend.instr"
    u);
#line 10949 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_4, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaelliptical_guide_gravity3 = coords_add(mcposaEndOfelement_4, mctc2);
  mctc1 = coords_sub(mcposaEndOfelement_4, mcposaelliptical_guide_gravity3);
  mcposrelliptical_guide_gravity3 = rot_apply(mcrotaelliptical_guide_gravity3, mctc1);
  mcDEBUG_COMPONENT("elliptical_guide_gravity3", mcposaelliptical_guide_gravity3, mcrotaelliptical_guide_gravity3)
  mccomp_posa[7] = mcposaelliptical_guide_gravity3;
  mccomp_posr[7] = mcposrelliptical_guide_gravity3;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component EndOfelement_3. */
  /* Setting parameters for component EndOfelement_3. */
  SIG_MESSAGE("EndOfelement_3 (Init:SetPar)");

  SIG_MESSAGE("EndOfelement_3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10969 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaelliptical_guide_gravity3, mcrotaEndOfelement_3);
  rot_transpose(mcrotaelliptical_guide_gravity3, mctr1);
  rot_mul(mcrotaEndOfelement_3, mctr1, mcrotrEndOfelement_3);
  mctc1 = coords_set(
#line 1378 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1378 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1378 "FullInstrument_v4.7_frontend.instr"
    length3 + 2 * u);
#line 10980 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaelliptical_guide_gravity3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEndOfelement_3 = coords_add(mcposaelliptical_guide_gravity3, mctc2);
  mctc1 = coords_sub(mcposaelliptical_guide_gravity3, mcposaEndOfelement_3);
  mcposrEndOfelement_3 = rot_apply(mcrotaEndOfelement_3, mctc1);
  mcDEBUG_COMPONENT("EndOfelement_3", mcposaEndOfelement_3, mcrotaEndOfelement_3)
  mccomp_posa[8] = mcposaEndOfelement_3;
  mccomp_posr[8] = mcposrEndOfelement_3;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component Curved_guide_2. */
  /* Setting parameters for component Curved_guide_2. */
  SIG_MESSAGE("Curved_guide_2 (Init:SetPar)");
#line 1382 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_w = startx2;
#line 1382 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_h = starty2;
#line 1383 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_r = curve_radius2;
#line 94 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Win = 0.04;
#line 1383 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_k = channels2;
#line 1383 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_d = 0.0001;
#line 1382 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_l = length2;
#line 1384 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_R0a = R02;
#line 1384 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Qca = Qc2;
#line 1384 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_alphaa = m2 * 0.24 + 2.1;
#line 1384 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_ma = m2;
#line 1384 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Wa = W2;
#line 1385 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_R0i = R02;
#line 1385 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Qci = Qc2;
#line 1385 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_alphai = m2 * 0.24 + 2.1;
#line 1385 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_mi = m2;
#line 1385 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Wi = W2;
#line 1386 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_R0s = R02;
#line 1386 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Qcs = Qc2;
#line 1386 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_alphas = m2 * 0.24 + 2.1;
#line 1386 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_ms = m2;
#line 1386 "FullInstrument_v4.7_frontend.instr"
  mccCurved_guide_2_Ws = W2;
#line 11038 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("Curved_guide_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1388 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1388 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1388 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11048 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_3, mcrotaCurved_guide_2);
  rot_transpose(mcrotaEndOfelement_3, mctr1);
  rot_mul(mcrotaCurved_guide_2, mctr1, mcrotrCurved_guide_2);
  mctc1 = coords_set(
#line 1387 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1387 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1387 "FullInstrument_v4.7_frontend.instr"
    u);
#line 11059 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_3, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaCurved_guide_2 = coords_add(mcposaEndOfelement_3, mctc2);
  mctc1 = coords_sub(mcposaEndOfelement_3, mcposaCurved_guide_2);
  mcposrCurved_guide_2 = rot_apply(mcrotaCurved_guide_2, mctc1);
  mcDEBUG_COMPONENT("Curved_guide_2", mcposaCurved_guide_2, mcrotaCurved_guide_2)
  mccomp_posa[9] = mcposaCurved_guide_2;
  mccomp_posr[9] = mcposrCurved_guide_2;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component EndOfelement_2. */
  /* Setting parameters for component EndOfelement_2. */
  SIG_MESSAGE("EndOfelement_2 (Init:SetPar)");

  SIG_MESSAGE("EndOfelement_2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1392 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1392 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1392 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11082 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaCurved_guide_2, mcrotaEndOfelement_2);
  rot_transpose(mcrotaCurved_guide_2, mctr1);
  rot_mul(mcrotaEndOfelement_2, mctr1, mcrotrEndOfelement_2);
  mctc1 = coords_set(
#line 1391 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1391 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1391 "FullInstrument_v4.7_frontend.instr"
    length2 + 2 * u);
#line 11093 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaCurved_guide_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEndOfelement_2 = coords_add(mcposaCurved_guide_2, mctc2);
  mctc1 = coords_sub(mcposaCurved_guide_2, mcposaEndOfelement_2);
  mcposrEndOfelement_2 = rot_apply(mcrotaEndOfelement_2, mctc1);
  mcDEBUG_COMPONENT("EndOfelement_2", mcposaEndOfelement_2, mcrotaEndOfelement_2)
  mccomp_posa[10] = mcposaEndOfelement_2;
  mccomp_posr[10] = mcposrEndOfelement_2;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component elliptical_guide_gravity1. */
  /* Setting parameters for component elliptical_guide_gravity1. */
  SIG_MESSAGE("elliptical_guide_gravity1 (Init:SetPar)");
#line 1396 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_l = length1;
#line 1397 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_xwidth = 2 * smallaxis_x1;
#line 1397 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_yheight = 2 * smallaxis_y1;
#line 1396 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_linxw = Linx1;
#line 1396 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_loutxw = Loutx1;
#line 1396 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_linyh = Liny1;
#line 1396 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_loutyh = Louty1;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_majorAxisxw = 0;
#line 165 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_minorAxisxw = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_majorAxisyh = 0;
#line 166 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_minorAxisyh = 0;
#line 167 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_majorAxisoffsetxw = 0;
#line 168 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_majorAxisoffsetyh = 0;
#line 1397 "FullInstrument_v4.7_frontend.instr"
  if("mid") strncpy(mccelliptical_guide_gravity1_dimensionsAt, "mid" ? "mid" : "", 16384); else mccelliptical_guide_gravity1_dimensionsAt[0]='\0';
#line 170 "FullInstrument_v4.7_frontend.instr"
  if("ellipse") strncpy(mccelliptical_guide_gravity1_option, "ellipse" ? "ellipse" : "", 16384); else mccelliptical_guide_gravity1_option[0]='\0';
#line 1398 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_R0 = R01;
#line 1398 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_Qc = Qc1;
#line 1398 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_alpha = m1 * 0.24 + 2.1;
#line 1398 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_m = m1;
#line 1398 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_W = W1;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_alpharight = -1;
#line 172 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_mright = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_alphaleft = -1;
#line 173 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_mleft = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_alphatop = -1;
#line 174 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_mtop = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_alphabottom = -1;
#line 175 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_mbottom = -1;
#line 176 "FullInstrument_v4.7_frontend.instr"
  if("on") strncpy(mccelliptical_guide_gravity1_verbose, "on" ? "on" : "", 16384); else mccelliptical_guide_gravity1_verbose[0]='\0';
#line 177 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_enableGravity = 1.0;
#line 178 "FullInstrument_v4.7_frontend.instr"
  mccelliptical_guide_gravity1_curvature = 0;
#line 11169 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("elliptical_guide_gravity1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1400 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1400 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1400 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11179 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_2, mcrotaelliptical_guide_gravity1);
  rot_transpose(mcrotaEndOfelement_2, mctr1);
  rot_mul(mcrotaelliptical_guide_gravity1, mctr1, mcrotrelliptical_guide_gravity1);
  mctc1 = coords_set(
#line 1399 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1399 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1399 "FullInstrument_v4.7_frontend.instr"
    u);
#line 11190 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaelliptical_guide_gravity1 = coords_add(mcposaEndOfelement_2, mctc2);
  mctc1 = coords_sub(mcposaEndOfelement_2, mcposaelliptical_guide_gravity1);
  mcposrelliptical_guide_gravity1 = rot_apply(mcrotaelliptical_guide_gravity1, mctc1);
  mcDEBUG_COMPONENT("elliptical_guide_gravity1", mcposaelliptical_guide_gravity1, mcrotaelliptical_guide_gravity1)
  mccomp_posa[11] = mcposaelliptical_guide_gravity1;
  mccomp_posr[11] = mcposrelliptical_guide_gravity1;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component EndOfelement_1. */
  /* Setting parameters for component EndOfelement_1. */
  SIG_MESSAGE("EndOfelement_1 (Init:SetPar)");

  SIG_MESSAGE("EndOfelement_1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11210 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaelliptical_guide_gravity1, mcrotaEndOfelement_1);
  rot_transpose(mcrotaelliptical_guide_gravity1, mctr1);
  rot_mul(mcrotaEndOfelement_1, mctr1, mcrotrEndOfelement_1);
  mctc1 = coords_set(
#line 1403 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1403 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1403 "FullInstrument_v4.7_frontend.instr"
    length1 + 2 * u);
#line 11221 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaelliptical_guide_gravity1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaEndOfelement_1 = coords_add(mcposaelliptical_guide_gravity1, mctc2);
  mctc1 = coords_sub(mcposaelliptical_guide_gravity1, mcposaEndOfelement_1);
  mcposrEndOfelement_1 = rot_apply(mcrotaEndOfelement_1, mctc1);
  mcDEBUG_COMPONENT("EndOfelement_1", mcposaEndOfelement_1, mcrotaEndOfelement_1)
  mccomp_posa[12] = mcposaEndOfelement_1;
  mccomp_posr[12] = mcposrEndOfelement_1;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component SlitGuideEnd. */
  /* Setting parameters for component SlitGuideEnd. */
  SIG_MESSAGE("SlitGuideEnd (Init:SetPar)");
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_xmin = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_xmax = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_ymin = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_ymax = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_radius = 0;
#line 1440 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_xwidth = mcipslitwidth;
#line 1440 "FullInstrument_v4.7_frontend.instr"
  mccSlitGuideEnd_yheight = slitheight;
#line 11249 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("SlitGuideEnd (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11256 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_1, mcrotaSlitGuideEnd);
  rot_transpose(mcrotaEndOfelement_1, mctr1);
  rot_mul(mcrotaSlitGuideEnd, mctr1, mcrotrSlitGuideEnd);
  mctc1 = coords_set(
#line 1441 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1441 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1441 "FullInstrument_v4.7_frontend.instr"
    0.001);
#line 11267 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSlitGuideEnd = coords_add(mcposaEndOfelement_1, mctc2);
  mctc1 = coords_sub(mcposaEndOfelement_1, mcposaSlitGuideEnd);
  mcposrSlitGuideEnd = rot_apply(mcrotaSlitGuideEnd, mctc1);
  mcDEBUG_COMPONENT("SlitGuideEnd", mcposaSlitGuideEnd, mcrotaSlitGuideEnd)
  mccomp_posa[13] = mcposaSlitGuideEnd;
  mccomp_posr[13] = mcposrSlitGuideEnd;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component SlitMonitor. */
  /* Setting parameters for component SlitMonitor. */
  SIG_MESSAGE("SlitMonitor (Init:SetPar)");
#line 1444 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_nx = 100;
#line 1445 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_ny = 100;
#line 1446 "FullInstrument_v4.7_frontend.instr"
  if("SlitMonitor.dat") strncpy(mccSlitMonitor_filename, "SlitMonitor.dat" ? "SlitMonitor.dat" : "", 16384); else mccSlitMonitor_filename[0]='\0';
#line 1447 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_xmin = -0.05;
#line 1448 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_xmax = 0.05;
#line 1449 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_ymin = -0.05;
#line 1450 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_ymax = 0.05;
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_xwidth = 0;
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_yheight = 0;
#line 1451 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_restore_neutron = 1;
#line 52 "FullInstrument_v4.7_frontend.instr"
  mccSlitMonitor_nowritefile = 0;
#line 11303 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("SlitMonitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11310 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaEndOfelement_1, mcrotaSlitMonitor);
  rot_transpose(mcrotaSlitGuideEnd, mctr1);
  rot_mul(mcrotaSlitMonitor, mctr1, mcrotrSlitMonitor);
  mctc1 = coords_set(
#line 1452 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1452 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1452 "FullInstrument_v4.7_frontend.instr"
    mcipvirtual_slit_d);
#line 11321 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaEndOfelement_1, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSlitMonitor = coords_add(mcposaEndOfelement_1, mctc2);
  mctc1 = coords_sub(mcposaSlitGuideEnd, mcposaSlitMonitor);
  mcposrSlitMonitor = rot_apply(mcrotaSlitMonitor, mctc1);
  mcDEBUG_COMPONENT("SlitMonitor", mcposaSlitMonitor, mcrotaSlitMonitor)
  mccomp_posa[14] = mcposaSlitMonitor;
  mccomp_posr[14] = mcposrSlitMonitor;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component MonoFLat. */
  /* Setting parameters for component MonoFLat. */
  SIG_MESSAGE("MonoFLat (Init:SetPar)");
#line 99 "FullInstrument_v4.7_frontend.instr"
  if("NULL") strncpy(mccMonoFLat_reflect, "NULL" ? "NULL" : "", 16384); else mccMonoFLat_reflect[0]='\0';
#line 99 "FullInstrument_v4.7_frontend.instr"
  if("NULL") strncpy(mccMonoFLat_transmit, "NULL" ? "NULL" : "", 16384); else mccMonoFLat_transmit[0]='\0';
#line 1468 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_zwidth = 0.05;
#line 1468 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_yheight = 0.025;
#line 1468 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_gap = 0.0005;
#line 1468 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_NH = 3;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_NV = 5;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_mosaich = 37;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_mosaicv = 37;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_r0 = 1;
#line 101 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_t0 = 1.0;
#line 101 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_Q = 1.8734;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_RV = RV_mono;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_RH = 0;
#line 1469 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_DM = 3.355;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_mosaic = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_width = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_height = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_verbose = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMonoFLat_order = 0;
#line 11375 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("MonoFLat (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1471 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1471 "FullInstrument_v4.7_frontend.instr"
    (thetaA1)*DEG2RAD,
#line 1471 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11385 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaSlitGuideEnd, mcrotaMonoFLat);
  rot_transpose(mcrotaSlitMonitor, mctr1);
  rot_mul(mcrotaMonoFLat, mctr1, mcrotrMonoFLat);
  mctc1 = coords_set(
#line 1470 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1470 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1470 "FullInstrument_v4.7_frontend.instr"
    1.6);
#line 11396 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaSlitGuideEnd, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMonoFLat = coords_add(mcposaSlitGuideEnd, mctc2);
  mctc1 = coords_sub(mcposaSlitMonitor, mcposaMonoFLat);
  mcposrMonoFLat = rot_apply(mcrotaMonoFLat, mctc1);
  mcDEBUG_COMPONENT("MonoFLat", mcposaMonoFLat, mcrotaMonoFLat)
  mccomp_posa[15] = mcposaMonoFLat;
  mccomp_posr[15] = mcposrMonoFLat;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component Mono. */
  /* Setting parameters for component Mono. */
  SIG_MESSAGE("Mono (Init:SetPar)");
#line 99 "FullInstrument_v4.7_frontend.instr"
  if("NULL") strncpy(mccMono_reflect, "NULL" ? "NULL" : "", 16384); else mccMono_reflect[0]='\0';
#line 99 "FullInstrument_v4.7_frontend.instr"
  if("NULL") strncpy(mccMono_transmit, "NULL" ? "NULL" : "", 16384); else mccMono_transmit[0]='\0';
#line 1475 "FullInstrument_v4.7_frontend.instr"
  mccMono_zwidth = 0.015;
#line 1475 "FullInstrument_v4.7_frontend.instr"
  mccMono_yheight = 0.015;
#line 1475 "FullInstrument_v4.7_frontend.instr"
  mccMono_gap = 0.002;
#line 1475 "FullInstrument_v4.7_frontend.instr"
  mccMono_NH = 17;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_NV = 11;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_mosaich = 42;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_mosaicv = 42;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_r0 = 1;
#line 101 "FullInstrument_v4.7_frontend.instr"
  mccMono_t0 = 1.0;
#line 101 "FullInstrument_v4.7_frontend.instr"
  mccMono_Q = 1.8734;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_RV = RV_mono;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_RH = RH_mono;
#line 1476 "FullInstrument_v4.7_frontend.instr"
  mccMono_DM = 3.355;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMono_mosaic = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMono_width = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMono_height = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMono_verbose = 0;
#line 102 "FullInstrument_v4.7_frontend.instr"
  mccMono_order = 0;
#line 11450 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("Mono (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1478 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1478 "FullInstrument_v4.7_frontend.instr"
    (thetaA1)*DEG2RAD,
#line 1478 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11460 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaSlitGuideEnd, mcrotaMono);
  rot_transpose(mcrotaMonoFLat, mctr1);
  rot_mul(mcrotaMono, mctr1, mcrotrMono);
  mctc1 = coords_set(
#line 1477 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1477 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1477 "FullInstrument_v4.7_frontend.instr"
    1.6);
#line 11471 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaSlitGuideEnd, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaMono = coords_add(mcposaSlitGuideEnd, mctc2);
  mctc1 = coords_sub(mcposaMonoFLat, mcposaMono);
  mcposrMono = rot_apply(mcrotaMono, mctc1);
  mcDEBUG_COMPONENT("Mono", mcposaMono, mcrotaMono)
  mccomp_posa[16] = mcposaMono;
  mccomp_posr[16] = mcposrMono;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component A2. */
  /* Setting parameters for component A2. */
  SIG_MESSAGE("A2 (Init:SetPar)");

  SIG_MESSAGE("A2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 1484 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD,
#line 1484 "FullInstrument_v4.7_frontend.instr"
    (thetaA2)*DEG2RAD,
#line 1484 "FullInstrument_v4.7_frontend.instr"
    (0)*DEG2RAD);
#line 11494 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaSlitGuideEnd, mcrotaA2);
  rot_transpose(mcrotaMono, mctr1);
  rot_mul(mcrotaA2, mctr1, mcrotrA2);
  mctc1 = coords_set(
#line 1483 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1483 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1483 "FullInstrument_v4.7_frontend.instr"
    1.6);
#line 11505 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaSlitGuideEnd, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaA2 = coords_add(mcposaSlitGuideEnd, mctc2);
  mctc1 = coords_sub(mcposaMono, mcposaA2);
  mcposrA2 = rot_apply(mcrotaA2, mctc1);
  mcDEBUG_COMPONENT("A2", mcposaA2, mcrotaA2)
  mccomp_posa[17] = mcposaA2;
  mccomp_posr[17] = mcposrA2;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component slit. */
  /* Setting parameters for component slit. */
  SIG_MESSAGE("slit (Init:SetPar)");
#line 1505 "FullInstrument_v4.7_frontend.instr"
  mccslit_xmin = -0.005;
#line 1506 "FullInstrument_v4.7_frontend.instr"
  mccslit_xmax = 0.005;
#line 1507 "FullInstrument_v4.7_frontend.instr"
  mccslit_ymin = -0.005;
#line 1508 "FullInstrument_v4.7_frontend.instr"
  mccslit_ymax = 0.005;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccslit_radius = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccslit_xwidth = 0;
#line 46 "FullInstrument_v4.7_frontend.instr"
  mccslit_yheight = 0;
#line 11533 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11540 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaA2, mcrotaslit);
  rot_transpose(mcrotaA2, mctr1);
  rot_mul(mcrotaslit, mctr1, mcrotrslit);
  mctc1 = coords_set(
#line 1509 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1509 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1509 "FullInstrument_v4.7_frontend.instr"
    1.5);
#line 11551 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaA2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit = coords_add(mcposaA2, mctc2);
  mctc1 = coords_sub(mcposaA2, mcposaslit);
  mcposrslit = rot_apply(mcrotaslit, mctc1);
  mcDEBUG_COMPONENT("slit", mcposaslit, mcrotaslit)
  mccomp_posa[18] = mcposaslit;
  mccomp_posr[18] = mcposrslit;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component SamplePos. */
  /* Setting parameters for component SamplePos. */
  SIG_MESSAGE("SamplePos (Init:SetPar)");

  SIG_MESSAGE("SamplePos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11571 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaA2, mcrotaSamplePos);
  rot_transpose(mcrotaslit, mctr1);
  rot_mul(mcrotaSamplePos, mctr1, mcrotrSamplePos);
  mctc1 = coords_set(
#line 1513 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1513 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1513 "FullInstrument_v4.7_frontend.instr"
    1.6);
#line 11582 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaA2, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSamplePos = coords_add(mcposaA2, mctc2);
  mctc1 = coords_sub(mcposaslit, mcposaSamplePos);
  mcposrSamplePos = rot_apply(mcrotaSamplePos, mctc1);
  mcDEBUG_COMPONENT("SamplePos", mcposaSamplePos, mcrotaSamplePos)
  mccomp_posa[19] = mcposaSamplePos;
  mccomp_posr[19] = mcposrSamplePos;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component psd_monitor_sample_pos. */
  /* Setting parameters for component psd_monitor_sample_pos. */
  SIG_MESSAGE("psd_monitor_sample_pos (Init:SetPar)");
#line 1516 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_nx = 101;
#line 1517 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_ny = 101;
#line 1518 "FullInstrument_v4.7_frontend.instr"
  if("psd_sampel_pos") strncpy(mccpsd_monitor_sample_pos_filename, "psd_sampel_pos" ? "psd_sampel_pos" : "", 16384); else mccpsd_monitor_sample_pos_filename[0]='\0';
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_xmin = -0.05;
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_xmax = 0.05;
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_ymin = -0.05;
#line 51 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_ymax = 0.05;
#line 1519 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_xwidth = 2 * mcipsampleRadius;
#line 1520 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_yheight = mcipsampleHeight;
#line 1521 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_restore_neutron = 1;
#line 52 "FullInstrument_v4.7_frontend.instr"
  mccpsd_monitor_sample_pos_nowritefile = 0;
#line 11618 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("psd_monitor_sample_pos (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11625 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotaSamplePos, mcrotapsd_monitor_sample_pos);
  rot_transpose(mcrotaSamplePos, mctr1);
  rot_mul(mcrotapsd_monitor_sample_pos, mctr1, mcrotrpsd_monitor_sample_pos);
  mctc1 = coords_set(
#line 1522 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1522 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1522 "FullInstrument_v4.7_frontend.instr"
    0);
#line 11636 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotaSamplePos, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_sample_pos = coords_add(mcposaSamplePos, mctc2);
  mctc1 = coords_sub(mcposaSamplePos, mcposapsd_monitor_sample_pos);
  mcposrpsd_monitor_sample_pos = rot_apply(mcrotapsd_monitor_sample_pos, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_sample_pos", mcposapsd_monitor_sample_pos, mcrotapsd_monitor_sample_pos)
  mccomp_posa[20] = mcposapsd_monitor_sample_pos;
  mccomp_posr[20] = mcposrpsd_monitor_sample_pos;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component e_monitor_sample_ps. */
  /* Setting parameters for component e_monitor_sample_ps. */
  SIG_MESSAGE("e_monitor_sample_ps (Init:SetPar)");
#line 1526 "FullInstrument_v4.7_frontend.instr"
  if("e_mon_sample") strncpy(mcce_monitor_sample_ps_filename, "e_mon_sample" ? "e_mon_sample" : "", 16384); else mcce_monitor_sample_ps_filename[0]='\0';
#line 53 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_xmin = -0.05;
#line 53 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_xmax = 0.05;
#line 53 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_ymin = -0.05;
#line 53 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_ymax = 0.05;
#line 1527 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_xwidth = 2 * mcipsampleRadius;
#line 1528 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_yheight = mcipsampleHeight;
#line 1529 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_Emin = 3.1;
#line 1530 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_Emax = 5.1;
#line 54 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_restore_neutron = 0;
#line 54 "FullInstrument_v4.7_frontend.instr"
  mcce_monitor_sample_ps_nowritefile = 0;
#line 11672 "./FullInstrument_v4.7_frontend.c"

  SIG_MESSAGE("e_monitor_sample_ps (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 11679 "./FullInstrument_v4.7_frontend.c"
  rot_mul(mctr1, mcrotapsd_monitor_sample_pos, mcrotae_monitor_sample_ps);
  rot_transpose(mcrotapsd_monitor_sample_pos, mctr1);
  rot_mul(mcrotae_monitor_sample_ps, mctr1, mcrotre_monitor_sample_ps);
  mctc1 = coords_set(
#line 1531 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1531 "FullInstrument_v4.7_frontend.instr"
    0,
#line 1531 "FullInstrument_v4.7_frontend.instr"
    0.001);
#line 11690 "./FullInstrument_v4.7_frontend.c"
  rot_transpose(mcrotapsd_monitor_sample_pos, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposae_monitor_sample_ps = coords_add(mcposapsd_monitor_sample_pos, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_sample_pos, mcposae_monitor_sample_ps);
  mcposre_monitor_sample_ps = rot_apply(mcrotae_monitor_sample_ps, mctc1);
  mcDEBUG_COMPONENT("e_monitor_sample_ps", mcposae_monitor_sample_ps, mcrotae_monitor_sample_ps)
  mccomp_posa[21] = mcposae_monitor_sample_ps;
  mccomp_posr[21] = mcposre_monitor_sample_ps;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
#define profile mccOrigin_profile
#define percent mccOrigin_percent
#define flag_save mccOrigin_flag_save
#define minutes mccOrigin_minutes
#line 57 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 11727 "./FullInstrument_v4.7_frontend.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component source. */
  SIG_MESSAGE("source (Init)");
#define mccompcurname  source
#define mccompcurtype  Source_gen4
#define mccompcurindex 2
#define flux_file mccsource_flux_file
#define xdiv_file mccsource_xdiv_file
#define ydiv_file mccsource_ydiv_file
#define p_in mccsource_p_in
#define lambda0 mccsource_lambda0
#define lambda02 mccsource_lambda02
#define L2P mccsource_L2P
#define lambda0b mccsource_lambda0b
#define lambda02b mccsource_lambda02b
#define L2Pb mccsource_L2Pb
#define lambda0c mccsource_lambda0c
#define lambda02c mccsource_lambda02c
#define L2Pc mccsource_L2Pc
#define pTable mccsource_pTable
#define pTable_x mccsource_pTable_x
#define pTable_y mccsource_pTable_y
#define radius mccsource_radius
#define dist mccsource_dist
#define xw mccsource_xw
#define yh mccsource_yh
#define E0 mccsource_E0
#define dE mccsource_dE
#define Lambda0 mccsource_Lambda0
#define dLambda mccsource_dLambda
#define I1 mccsource_I1
#define h mccsource_h
#define w mccsource_w
#define gaussian mccsource_gaussian
#define verbose mccsource_verbose
#define T1 mccsource_T1
#define flux_file_perAA mccsource_flux_file_perAA
#define flux_file_log mccsource_flux_file_log
#define Lmin mccsource_Lmin
#define Lmax mccsource_Lmax
#define Emin mccsource_Emin
#define Emax mccsource_Emax
#define T2 mccsource_T2
#define I2 mccsource_I2
#define T3 mccsource_T3
#define I3 mccsource_I3
#define length mccsource_length
#define phi_init mccsource_phi_init
#define theta_init mccsource_theta_init
#define HEtailA mccsource_HEtailA
#define HEtailL0 mccsource_HEtailL0
#line 190 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
{
  double source_area, k;

  /* spectrum characteristics */
  if (flux_file && strlen(flux_file) > 0) {
    if (Table_Read(&pTable, flux_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read file %s\n", NAME_CURRENT_COMP, flux_file));
    /* put table in Log scale */
    int i;
    if (pTable.columns < 2) exit(fprintf(stderr, "Source_gen: %s: Flux file %s should contain at least 2 columns\n", NAME_CURRENT_COMP, flux_file));
    double table_lmin=FLT_MAX, table_lmax=-FLT_MAX;
    double tmin=FLT_MAX, tmax=-FLT_MAX;
    for (i=0; i<pTable.rows; i++) {
      double val = Table_Index(pTable, i,1);
      val = Table_Index(pTable, i,0); /* lambda */
      if (val > tmax) tmax=val;
      if (val < tmin) tmin=val;
    }
    for (i=0; i<pTable.rows; i++) {
      double val = Table_Index(pTable, i,1);
      if (val < 0) fprintf(stderr, "Source_gen: %s: File %s has negative flux at row %i\n", NAME_CURRENT_COMP, flux_file, i+1);
      if (flux_file_log)
        val = log(val > 0 ? val : tmin/10);
      Table_SetElement(&pTable, i, 1, val);
      val = Table_Index(pTable, i,0); /* lambda */
      if (val > table_lmax) table_lmax=val;
      if (val < table_lmin) table_lmin=val;
    }
    if (!Lmin && !Lmax && !Lambda0 && !dLambda && !E0 && !dE && !Emin && !Emax) {
      Lmin = table_lmin; Lmax = table_lmax;
    }
    if (Lmax > table_lmax) {
      if (verbose) fprintf(stderr, "Source_gen: %s: Maximum wavelength %g is beyond table range upper limit %g. Constraining.\n", NAME_CURRENT_COMP, Lmax, table_lmax);
      Lmax = table_lmax;
    }
    if (Lmin < table_lmin) {
      if (verbose) fprintf(stderr, "Source_gen: %s: Minimum wavelength %g is below table range lower limit %g. Constraining.\n", NAME_CURRENT_COMP, Lmin, table_lmin);
      Lmin = table_lmin;
    }
  } else
  {
    k  = 1.38066e-23; /* k_B */
    if (T1 > 0)
    {
      lambda0  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T1);
      lambda02 = lambda0*lambda0;
      L2P      = 2*lambda02*lambda02;
    }
    else
      { lambda0 = Lambda0; }

    if (T2 > 0)
    {
      lambda0b  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T2);
      lambda02b = lambda0b*lambda0b;
      L2Pb      = 2*lambda02b*lambda02b;
    }
    else
      { lambda0b = Lambda0; }

    if (T3 > 0)
    {
      lambda0c  = 1.0e10*sqrt(HBAR*HBAR*4.0*PI*PI/2.0/MNEUTRON/k/T3);
      lambda02c = lambda0c*lambda0c;
      L2Pc      = 2*lambda02c*lambda02c;
    }
    else
      { lambda0c = Lambda0; }
  }

  /* now read position-divergence files, if any */
  if (xdiv_file && strlen(xdiv_file) > 0) {
    int i,j;
    if (Table_Read(&pTable_x, xdiv_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read file %s\n", NAME_CURRENT_COMP, xdiv_file));
    pTable_xsum = 0;
    for (i=0; i<pTable_x.rows; i++)
      for (j=0; j<pTable_x.columns; j++)
        pTable_xsum += Table_Index(pTable_x, i,j);

    /* now extract limits */
    char **parsing;
    char xylimits[1024];
    strcpy(xylimits, "");
    parsing = Table_ParseHeader(pTable_x.header,
      "xlimits", "xylimits",
      NULL);

    if (parsing) {
      if (parsing[0])  strcpy(xylimits, str_dup_numeric(parsing[0]));
      if (parsing[1] && !strlen(xylimits))
                       strcpy(xylimits, str_dup_numeric(parsing[1]));
      for (i=0; i<=1; i++) {
        if (parsing[i]) free(parsing[i]);
      }
      free(parsing);
    }
    i = sscanf(xylimits, "%lg %lg %lg %lg",
      &(pTable_xmin),  &(pTable_xmax),
      &(pTable_dxmin), &(pTable_dxmax));
    if (i != 2 && i != 4 && verbose)
      fprintf(stderr, "Source_gen: %s: invalid xylimits '%s' from file %s. extracted %i values\n",
        NAME_CURRENT_COMP, xylimits, xdiv_file, i);

    if (!w) w=pTable_xmax-pTable_xmin;
    if (!xw && !dist) xw=fabs(pTable_dxmax-pTable_dxmin);
  }

  if (ydiv_file && strlen(ydiv_file) > 0) {
    int i,j;
    if (Table_Read(&pTable_y, ydiv_file, 1) <= 0) /* read 1st block data from file into pTable */
      exit(fprintf(stderr, "Source_gen: %s: can not read file %s\n", NAME_CURRENT_COMP, ydiv_file));
    pTable_ysum = 0;
    for (i=0; i<pTable_y.rows; i++)
      for (j=0; j<pTable_y.columns; j++)
        pTable_ysum += Table_Index(pTable_y, i,j);

    /* now extract limits */
    char **parsing;
    char xylimits[1024];
    strcpy(xylimits, "");
    parsing = Table_ParseHeader(pTable_y.header,
      "xlimits", "xylimits",
      NULL);

    if (parsing) {
      if (parsing[0])  strcpy(xylimits,str_dup_numeric(parsing[0]));
      if (parsing[1] && !strlen(xylimits))
                       strcpy(xylimits,str_dup_numeric(parsing[1]));
      for (i=0; i<=1; i++) {
        if (parsing[i]) free(parsing[i]);
      }
      free(parsing);
    }
    i = sscanf(xylimits, "%lg %lg %lg %lg",
      &(pTable_ymin),  &(pTable_ymax),
      &(pTable_dymin), &(pTable_dymax));
    if (i != 2 && i != 4 && verbose)
      fprintf(stderr, "Source_gen: %s: invalid xylimits '%s' from file %s. extracted %i values\n",
        NAME_CURRENT_COMP, xylimits, ydiv_file, i);
    if (!h)  h=pTable_ymax-pTable_ymin;
    if (!yh && !dist) yh=fabs(pTable_dymax-pTable_dymin);
  }

  /* tests for parameter values */
  if (Emin < 0 || Emax < 0 || Lmin < 0 || Lmax < 0 || E0 < 0 || dE < 0 || Lambda0 < 0 || dLambda < 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Negative average\n"
                   "            or range values for wavelength or energy encountered\n",
                   NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((Emin == 0 && Emax > 0) || (dE > 0 && dE >= E0))
  {
    fprintf(stderr,"Source_gen: %s: Error: minimal energy cannot be less or equal zero\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((Emax >= Emin) && (Emin > 0))
  { E0 = (Emax+Emin)/2;
    dE = (Emax-Emin)/2;
  }
  if ((E0 > dE) && (dE >= 0))
  {
    Lmin = sqrt(81.81/(E0+dE)); /* Angstroem */
    Lmax = sqrt(81.81/(E0-dE));
  }
  if (Lmax > 0)
  { Lambda0 = (Lmax+Lmin)/2;
    dLambda = (Lmax-Lmin)/2;
  }
  if ((Lambda0 < dLambda) || (dLambda < 0))
  { fprintf(stderr,"Source_gen: %s: Error: Wavelength range %.3f +/- %.3f AA calculated \n",
      NAME_CURRENT_COMP, Lambda0, dLambda);
    fprintf(stderr,"- whole wavelength range must be >= 0 \n");
    fprintf(stderr,"- range must be > 0; otherwise intensity gets zero, use other sources in this case \n\n");
    exit(-1);
  }

  radius = fabs(radius); w=fabs(w); h=fabs(h);  I1=fabs(I1);
  Lambda0=fabs(Lambda0); dLambda=fabs(dLambda);
  xw = fabs(xw); yh=fabs(yh); dist=fabs(dist);

  if (dist == 0)
  {
    fprintf(stderr,"Source_gen: %s: warning: focusing distance is null.\n"
                   "            xw and yh interpreted as full divergence in [deg]\n",
                   NAME_CURRENT_COMP);
  }
  Lmin = Lambda0 - dLambda; /* Angstroem */
  Lmax = Lambda0 + dLambda;

  /* compute initial weight factor p_in to get [n/s] */
  if ((I1 > 0  && T1 >= 0) || (flux_file && strlen(flux_file) > 0))
  { /* the I1,2,3 are usually in [n/s/cm2/st/AA] */
    if (radius)
      source_area = radius*radius*PI*1e4; /* circular cm^2 */
    else
      source_area = h*w*1e4; /* square cm^2 */
    p_in  = source_area; /* cm2 */
    p_in *= (Lmax-Lmin); /* AA. 1 bin=AA/n */
    if (flux_file && strlen(flux_file) && !flux_file_perAA)  p_in *= pTable.rows/(Lmax-Lmin);
  }
  else
    p_in = (I1 > 0? I1 : 1)/4/PI; /* Small angle approx. */
  p_in /= mcget_ncount();
  if (!T1 && I1) p_in *= I1;

  if (radius == 0 && h == 0 && w == 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Please specify source geometry (radius, h, w)\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }
  if (xw*yh == 0)
  {
    fprintf(stderr,"Source_gen: %s: Error: Please specify source target (xw, yh)\n",
      NAME_CURRENT_COMP);
    exit(-1);
  }

  if (verbose)
  {
    printf("Source_gen: component %s ", NAME_CURRENT_COMP);
    if ((h == 0) || (w == 0))
      printf("(disk, radius=%g)", radius);
    else
      printf("(square %g x %g)",h,w);
    printf("\n            spectra ");
    printf("%.3f to %.3f AA (%.3f to %.3f meV)", Lmin, Lmax, 81.81/Lmax/Lmax, 81.81/Lmin/Lmin);
    if (gaussian)
      printf(", gaussian divergence beam");
    printf("\n");
    if (flux_file && strlen(flux_file) > 0)
    { printf("  File %s for flux distribution used. Flux is dPhi/dLambda in [n/s/AA]. \n", flux_file);
      Table_Info(pTable);
    }
    else if (T1>=0 && I1)
    { if (T1 != 0)
        printf("            T1=%.1f K (%.3f AA)", T1, lambda0);
      if (T2*I2 != 0)
        printf(", T2=%.1f K (%.3f AA)", T2, lambda0b);
      if (T3*I3 != 0)
        printf(", T3=%.1f K (%.3f AA)", T3, lambda0c);
      if (T1) printf("\n");
      printf("  Flux is dPhi/dLambda in [n/s/cm2].\n");
    }
    else
    { printf("  Flux is Phi in [n/s].\n");
    }
    if (xdiv_file && strlen(xdiv_file) > 0)
      printf("  File %s x=[%g:%g] [m] xdiv=[%g:%g] [deg] used as horizontal phase space distribution.\n", xdiv_file, pTable_xmin, pTable_xmax, pTable_dxmin, pTable_dxmax);
    if (ydiv_file && strlen(ydiv_file) > 0)
      printf("  File %s y=[%g:%g] [m] ydiv=[%g:%g] [deg] used as vertical phase space distribution.\n", ydiv_file, pTable_ymin, pTable_ymax, pTable_dymin, pTable_dymax);
  }
  else
    if (verbose == -1)
      printf("Source_gen: component %s unactivated", NAME_CURRENT_COMP);
}
#line 12050 "./FullInstrument_v4.7_frontend.c"
#undef HEtailL0
#undef HEtailA
#undef theta_init
#undef phi_init
#undef length
#undef I3
#undef T3
#undef I2
#undef T2
#undef Emax
#undef Emin
#undef Lmax
#undef Lmin
#undef flux_file_log
#undef flux_file_perAA
#undef T1
#undef verbose
#undef gaussian
#undef w
#undef h
#undef I1
#undef dLambda
#undef Lambda0
#undef dE
#undef E0
#undef yh
#undef xw
#undef dist
#undef radius
#undef pTable_y
#undef pTable_x
#undef pTable
#undef L2Pc
#undef lambda02c
#undef lambda0c
#undef L2Pb
#undef lambda02b
#undef lambda0b
#undef L2P
#undef lambda02
#undef lambda0
#undef p_in
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component StartOfGuide. */
  SIG_MESSAGE("StartOfGuide (Init)");

  /* Initializations for component elliptical_guide_gravity5. */
  SIG_MESSAGE("elliptical_guide_gravity5 (Init)");
#define mccompcurname  elliptical_guide_gravity5
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 4
#define mvaluesright mccelliptical_guide_gravity5_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity5_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity5_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity5_mvaluesbottom
#define seglength mccelliptical_guide_gravity5_seglength
#define guideInfo mccelliptical_guide_gravity5_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity5_latestParticleCollision
#define Gx mccelliptical_guide_gravity5_Gx
#define Gy mccelliptical_guide_gravity5_Gy
#define Gz mccelliptical_guide_gravity5_Gz
#define Gx0 mccelliptical_guide_gravity5_Gx0
#define Gy0 mccelliptical_guide_gravity5_Gy0
#define Gz0 mccelliptical_guide_gravity5_Gz0
#define Circ mccelliptical_guide_gravity5_Circ
#define dynamicalSegLength mccelliptical_guide_gravity5_dynamicalSegLength
#define l mccelliptical_guide_gravity5_l
#define xwidth mccelliptical_guide_gravity5_xwidth
#define yheight mccelliptical_guide_gravity5_yheight
#define linxw mccelliptical_guide_gravity5_linxw
#define loutxw mccelliptical_guide_gravity5_loutxw
#define linyh mccelliptical_guide_gravity5_linyh
#define loutyh mccelliptical_guide_gravity5_loutyh
#define majorAxisxw mccelliptical_guide_gravity5_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity5_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity5_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity5_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity5_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity5_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity5_dimensionsAt
#define option mccelliptical_guide_gravity5_option
#define R0 mccelliptical_guide_gravity5_R0
#define Qc mccelliptical_guide_gravity5_Qc
#define alpha mccelliptical_guide_gravity5_alpha
#define m mccelliptical_guide_gravity5_m
#define W mccelliptical_guide_gravity5_W
#define alpharight mccelliptical_guide_gravity5_alpharight
#define mright mccelliptical_guide_gravity5_mright
#define alphaleft mccelliptical_guide_gravity5_alphaleft
#define mleft mccelliptical_guide_gravity5_mleft
#define alphatop mccelliptical_guide_gravity5_alphatop
#define mtop mccelliptical_guide_gravity5_mtop
#define alphabottom mccelliptical_guide_gravity5_alphabottom
#define mbottom mccelliptical_guide_gravity5_mbottom
#define verbose mccelliptical_guide_gravity5_verbose
#define enableGravity mccelliptical_guide_gravity5_enableGravity
#define curvature mccelliptical_guide_gravity5_curvature
#line 817 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
	///////////////////////////////////////////////////////////////////////////
	/////////////// Test user input
	///////////////////////////////////////////////////////////////////////////

	if(strcmp(verbose,"on") == 0)
		guideInfo.verboseSetting = 1;
	else guideInfo.verboseSetting = 0;
	

	guideInfo.R0Arr[RightSide] 		= R0;
	guideInfo.QcArr[RightSide] 		= Qc;
	guideInfo.alphaArr[RightSide] 	= alpharight;
	guideInfo.mArr[RightSide] 		= mright;
	guideInfo.WArr[RightSide] 		= W; 

	guideInfo.R0Arr[TopSide] 	 = R0; 
	guideInfo.QcArr[TopSide] 	 = Qc; 
	guideInfo.alphaArr[TopSide]  = alphatop; 
	guideInfo.mArr[TopSide] 	 = mtop;
	guideInfo.WArr[TopSide] 	 = W; 

	guideInfo.R0Arr[LeftSide] 		= R0; 
	guideInfo.QcArr[LeftSide] 		= Qc; 
	guideInfo.alphaArr[LeftSide] 	= alphaleft; 
	guideInfo.mArr[LeftSide] 		= mleft; 
	guideInfo.WArr[LeftSide] 		= W; 

	guideInfo.R0Arr[BottomSide] 	= R0; 
	guideInfo.QcArr[BottomSide] 	= Qc; 
	guideInfo.alphaArr[BottomSide] 	= alphabottom; 
	guideInfo.mArr[BottomSide] 		= mbottom; 
	guideInfo.WArr[BottomSide] 		= W; 

	int sides;
	for (sides = RightSide; sides <= BottomSide; sides++){
		if (guideInfo.alphaArr[sides] == -1) guideInfo.alphaArr[sides] = alpha;
		if (guideInfo.mArr[sides] == -1)	 guideInfo.mArr[sides] = m;
	}
	
	// Test user input for illegal values
	int inputErrors = 0;
	// Lower or equal to zero
	if(l <= 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("length",guideInfo.verboseSetting);
	if(guideInfo.alphaArr[TopSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphatop",guideInfo.verboseSetting);
	if(guideInfo.mArr[TopSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mtop",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[BottomSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphabottom",guideInfo.verboseSetting);
	if(guideInfo.mArr[BottomSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mbottom",guideInfo.verboseSetting);

	if(guideInfo.alphaArr[RightSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alpharight",guideInfo.verboseSetting);
	if(guideInfo.mArr[RightSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mright",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[LeftSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphaleft",guideInfo.verboseSetting);
	if(guideInfo.mArr[LeftSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mleft",guideInfo.verboseSetting);

	// Focal points result in hyperbola instead of an ellipse
	if(l <= -linxw-loutxw) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linw","loutw",linxw,loutxw,guideInfo.verboseSetting);
	if(l <= -linyh-loutyh) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linh","louth",linyh,loutyh,guideInfo.verboseSetting);

	if( strcmp(dimensionsAt,"entrance") != 0 
	 && strcmp(dimensionsAt,"mid") != 0 
	 && strcmp(dimensionsAt,"exit") != 0){
	 	inputErrors += 1;
		printf("dimensionsAt were given an incorrect input."
			   "Input must be string containing \"entrance\",\"mid\" or \"exit\" \n");
	}


	// Terminate program if any input errors occurred
	if(inputErrors != 0 ){
		exit(printf("\nCRITICAL ERROR(S) IN COMPONENT %s" 
				" CONSIDER CHECKING USER INPUT AS %d INPUT ERRORS WAS FOUND.\n",
				NAME_CURRENT_COMP,inputErrors) );
	}


	///////////////////////////////////////////////////////////////////////////
	/////////////// Calculate intern guide values from user input
	///////////////////////////////////////////////////////////////////////////

	/*  Calculate the foci line for the ellipses.
		These can be used to calculate the axes of the ellipses
		using pyth and defination of the ellipse that says distance 
		between the foci and every point on the ellipse is constant.
	*/
	int directDefination = 0;

	if( majorAxisyh != 0 || minorAxisyh != 0 
	 || majorAxisxw != 0 || minorAxisxw != 0)
	{
		directDefination = 1;
		guideInfo.Length = l;	
		
		guideInfo.ellipseMajorAxis[RightSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[RightSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[RightSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[TopSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[TopSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[TopSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[LeftSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[LeftSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[LeftSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[BottomSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[BottomSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[BottomSide] = 0;

		guideInfo.entranceHorizontalWidth =
			2*sqrt(1 - (majorAxisoffsetyh*majorAxisoffsetyh)
			/(majorAxisyh*majorAxisyh) )*minorAxisyh;
		guideInfo.entranceVerticalWidth = 	
			2*sqrt(1 - (majorAxisoffsetxw*majorAxisoffsetxw)
			/(majorAxisxw*majorAxisxw) )*minorAxisxw;
	}

	if ( strcmp(option,"ellipse") == 0 && directDefination == 0)
	{
		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( linxw*linxw + xwidth*xwidth*0.25)
				+ sqrt( (l + loutxw)*(l + loutxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( linyh*linyh + yheight*yheight*0.25)
				+ sqrt( (l + loutyh)*(l + loutyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"exit") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( loutxw*loutxw + xwidth*xwidth*0.25)
				+ sqrt( (l + linxw)*(l + linxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( loutyh*loutyh + yheight*yheight*0.25)
				+ sqrt( (l + linyh)*(l + linyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] =(l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"mid") == 0 ){ 

			guideInfo.Length = l;
			
			guideInfo.ellipseMajorAxis[RightSide] = 
				sqrt( (linxw+l+loutxw)*(linxw+l+loutxw)/4+xwidth*xwidth/4);
			guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = 
				sqrt( (linyh+l+loutyh)*(linyh+l+loutyh)/4+yheight*yheight/4);
			guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];

		}
	}

	guideInfo.entranceHorizontalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[RightSide]
			*guideInfo.ellipseMajorOffset[RightSide]
			/(guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] ) )
		*guideInfo.ellipseMinorAxis[RightSide];
	guideInfo.entranceVerticalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[TopSide]
			*guideInfo.ellipseMajorOffset[TopSide]
			/(guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] ) )
		*guideInfo.ellipseMinorAxis[TopSide];


	if ( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		exit( printf("Critical error in %s; the option for option = halfellipse is currently disabled.",NAME_CURRENT_COMP) );

		double used_focal_vertical;
		double used_focal_horizontal;
		double major_offset_horizontal = 0;
		double major_offset_vertical = 0;

		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+linyh)*(l+linyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+linxw)*(l+linxw) ); 
			major_offset_vertical = l;	
			major_offset_horizontal = l;
		}
		else {	
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+loutyh)*(l+loutyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+loutxw)*(l+loutxw) ); 
		}

		guideInfo.Length = l;

		guideInfo.ellipseMajorAxis[RightSide] = used_focal_horizontal;
		guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

		guideInfo.ellipseMajorOffset[RightSide] = major_offset_horizontal;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = 
								guideInfo.ellipseMajorAxis[RightSide];
		guideInfo.ellipseMinorAxis[LeftSide] = 
								guideInfo.ellipseMinorAxis[RightSide];
		guideInfo.ellipseMajorOffset[LeftSide] = 
								guideInfo.ellipseMajorOffset[RightSide];
		guideInfo.ellipseMinorOffset[LeftSide] = 
								guideInfo.ellipseMinorOffset[RightSide];

		guideInfo.ellipseMajorAxis[TopSide] = used_focal_vertical;
		guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

		guideInfo.ellipseMajorOffset[TopSide] = major_offset_vertical;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = 
								guideInfo.ellipseMajorAxis[TopSide];
		guideInfo.ellipseMinorAxis[BottomSide] = 
								guideInfo.ellipseMinorAxis[TopSide];
		guideInfo.ellipseMajorOffset[BottomSide] = 
								guideInfo.ellipseMajorOffset[TopSide];
		guideInfo.ellipseMinorOffset[BottomSide] = 
								guideInfo.ellipseMinorOffset[TopSide];
	}

	// Applies the properties of the mirrors in the guide given by the user.
	// These variables are used in the reflection functions.

	
	// Sets the mirror type of the guides mirrors
	// These variables are used in the collision functions 
	// to find the type of collision

	// guideInfo.OuterSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[BottomSide] = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_bottom_side_mirror,guideInfo.verboseSetting);

	// guideInfo.InnerSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[BottomSide] =	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_bottom_side_mirror,guideInfo.verboseSetting);

	// Give a warning if all side of the guide is turned off, 
	// as the guide is essentially turned off
	if(    guideInfo.OuterSide[RightSide] 	== 1 
		&& guideInfo.OuterSide[TopSide] 	== 1
		&& guideInfo.OuterSide[LeftSide] 	== 1 
		&& guideInfo.OuterSide[BottomSide]	== 1
		&& guideInfo.InnerSide[RightSide] 	== 1 
		&& guideInfo.InnerSide[TopSide] 	== 1
		&& guideInfo.InnerSide[LeftSide] 	== 1 
		&& guideInfo.InnerSide[BottomSide]  == 1)
		printf("Warning: In %s all the sides of the guide has been disabled,"
				" so it not possible for any particle"
				" to collide with the guide, consider"
				" disabling this component",NAME_CURRENT_COMP);
	
	if(guideInfo.mArr[RightSide] <= 0)  guideInfo.InnerSide[RightSide] = 
														MirrorTypeabsorption;
	if(guideInfo.mArr[TopSide] <= 0) 	 guideInfo.InnerSide[TopSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[LeftSide] <= 0) 	 guideInfo.InnerSide[LeftSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[BottomSide] <= 0) guideInfo.InnerSide[BottomSide] = 
														MirrorTypeabsorption;
	/* if(directDefination == 0){ */
	/* 	guideInfo.entranceHorizontalWidth = xwidth; */
	/* 	guideInfo.entranceVerticalWidth = yheight; */
	/* } */
	
	if( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		guideInfo.entranceHorizontalWidth =
			(guideInfo.ellipseMinorAxis[RightSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[RightSide]
						*guideInfo.ellipseMajorOffset[RightSide] )
			/( guideInfo.ellipseMajorAxis[RightSide]
			  *guideInfo.ellipseMajorAxis[RightSide] ) )
			+ guideInfo.ellipseMinorOffset[RightSide] )*2;
		guideInfo.entranceVerticalWidth = 
			(guideInfo.ellipseMinorAxis[TopSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[TopSide]
						*guideInfo.ellipseMajorOffset[TopSide] )
			/( guideInfo.ellipseMajorAxis[TopSide]
			  *guideInfo.ellipseMajorAxis[TopSide] ) )
			+ guideInfo.ellipseMinorOffset[TopSide] )*2;
	}


	guideInfo.EnclosingBoxOn = 0;

	/*
	double DefaultArray1[8] = { 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0};
	double DefaultArray2[8] = { 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0};
	double DefaultArray3[8] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0}; 

	guideInfo.EnclosingBoxOn = 0;
	double *xinput;
	if ( xInput != NULL ){ xinput = xInput; guideInfo.EnclosingBoxOn = 1; }
	else { xinput = DefaultArray1; }
	double *yinput;
	if ( yInput != NULL ){ yinput = yInput; guideInfo.EnclosingBoxOn = 1;}
	else { yinput = DefaultArray2; }
	double *zinput;
	if ( zInput != NULL ){ zinput = zInput; guideInfo.EnclosingBoxOn = 1;}
	else { zinput = DefaultArray3; }
	*/

	/*
	double xarray[8] ={ guideInfo.ellipseMinorAxis[0]*xinput[0],
						guideInfo.ellipseMinorAxis[2]*xinput[1],
						guideInfo.ellipseMinorAxis[2]*xinput[2],
						guideInfo.ellipseMinorAxis[0]*xinput[3],
						guideInfo.ellipseMinorAxis[0]*xinput[4],
						guideInfo.ellipseMinorAxis[2]*xinput[5],
						guideInfo.ellipseMinorAxis[2]*xinput[6],
						guideInfo.ellipseMinorAxis[0]*xinput[7] };
	double yarray[8] ={ guideInfo.ellipseMinorAxis[1]*yinput[0],
						guideInfo.ellipseMinorAxis[1]*yinput[1],
						guideInfo.ellipseMinorAxis[3]*yinput[2],
						guideInfo.ellipseMinorAxis[3]*yinput[3],
						guideInfo.ellipseMinorAxis[1]*yinput[4],
						guideInfo.ellipseMinorAxis[1]*yinput[5],
						guideInfo.ellipseMinorAxis[3]*yinput[6],
						guideInfo.ellipseMinorAxis[3]*yinput[7] };
	double zarray[8] ={ guideInfo.Length/2*zinput[0]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[1]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[2]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[3]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[4]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[5]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[6]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[7]+guideInfo.Length/2 };
	int i = 0;
	for(i = 0; i < 8; i++){
		guideInfo.xArray[i] = xarray[i];
		guideInfo.yArray[i] = yarray[i];
		guideInfo.zArray[i] = zarray[i];
	}
	*/

	guideInfo.exitVerticalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide]) 
					)/(guideInfo.ellipseMajorAxis[BottomSide] 
					*guideInfo.ellipseMajorAxis[BottomSide] )
		)*guideInfo.ellipseMinorAxis[BottomSide];

	guideInfo.exitHorizontalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide]) 
					)/(guideInfo.ellipseMajorAxis[RightSide] 
					*guideInfo.ellipseMajorAxis[RightSide] ) 
		)*guideInfo.ellipseMinorAxis[RightSide];

    //////////////////segmentation of m values

	// Are the arrays empty?
	if(mvaluesright != NULL || mvaluesleft != NULL 
	|| mvaluestop != NULL || mvaluesbottom != NULL)
	{
		guideInfo.enableSegments = 1;

		guideInfo.numberOfSegments = sizeof(mvaluesright)/sizeof(mvaluesright[0]);

		//printf("Length is %i\n",guideInfo.numberOfSegments);
		guideInfo.mValuesright = mvaluesright;
		guideInfo.mValuesleft = mvaluesleft;
		guideInfo.mValuestop = mvaluestop;
		guideInfo.mValuesbottom = mvaluesbottom;
		//printf("Seglength ... %f %f %f\n",seglength[0],seglength[1],seglength[2]);

		// Are the arrays of equal length?
		if(seglength == NULL){
			dynamicalSegLength = 
					realloc(dynamicalSegLength,
							guideInfo.numberOfSegments*sizeof(double) );
			int i;
			for (i = 0; i < guideInfo.numberOfSegments; ++i){
				dynamicalSegLength[i] = 
						guideInfo.Length/guideInfo.numberOfSegments;
				}
			guideInfo.segLength = dynamicalSegLength;
		}
		else guideInfo.segLength = seglength;

		if( guideInfo.numberOfSegments != sizeof(mvaluesright)/sizeof(mvaluesright[0]) 
		 || guideInfo.numberOfSegments != sizeof(mvaluesleft)/sizeof(mvaluesleft[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluestop)/sizeof(mvaluestop[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0])
		 || (guideInfo.segLength == NULL
		  & guideInfo.numberOfSegments != sizeof(seglength)/sizeof(guideInfo.segLength[0]) 
		  ) ) {
			
			printf("Error in userinput inside %s, the length of the arrays"
				   " mvalues and seglength are not equal\n",NAME_CURRENT_COMP);
			printf("The length of the arrays are: mValuesright is %lu,"
				   " mvaluesleft is %lu, mvaluestop is %lu, mvaluesbottom is"
				   " %lu and seglength is %lu and should be %d \n; Above assume that the arrays are using double \n",
					sizeof(mvaluesright)/sizeof(double),
					sizeof(mvaluesleft)/sizeof(double),
					sizeof(mvaluestop)/sizeof(double),
					sizeof(mvaluesbottom)/sizeof(double),
					sizeof(guideInfo.segLength)/sizeof(double),
					guideInfo.numberOfSegments	
					);
				
			if ( guideInfo.verboseSetting )	{
				int i;
				
				printf("The Values of mvaluesright is: [");
				for(i=0; i < sizeof(mvaluesright)/sizeof(mvaluesright[0]); i++) 	
						printf("%e,",guideInfo.mValuesright[i] );
				printf("]\n");
				
				printf("The Values of mvaluesleft is: [");
				for(i=0; i < sizeof(mvaluesleft)/sizeof(mvaluesleft[0]); i++) 	
						printf("%e,",guideInfo.mValuesleft[i] );
				printf("]\n");
				
				printf("The Values of mvaluestop is: [");
				for(i=0; i < sizeof(mvaluestop)/sizeof(mvaluestop[0]); i++) 	
						printf("%e,",guideInfo.mValuestop[i] );
				printf("]\n");

				printf("The Values of mvaluesbottom is: [");
				for(i=0; i < sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0]); i++) 
						printf("%e,",guideInfo.mValuesbottom[i] );
				printf("]\n");

				printf("The Values of seglength is: [");
				for(i=0; i < sizeof(guideInfo.segLength)/sizeof(guideInfo.segLength[0]); i++) 
						printf("%e,",guideInfo.segLength[i]);
				printf("]\n");
			}
			exit( printf("Exit due to critical error in userinput for the"
						" component %s, consider having a look at the input"
						" for following: mvaluesright,mvaluesleft,mvaluestop,"
						"mvaluesbottom and/or seglength.",NAME_CURRENT_COMP) );
		}
		//
		double sumOfelements=0;
		int i;
		for(i=0;i< guideInfo.numberOfSegments; i++) {
			sumOfelements += guideInfo.segLength[i];
		}
		if ( guideInfo.verboseSetting 
		&& fabs(sumOfelements-guideInfo.Length) > 1e-9 )
			printf("Error in userinput inside %s, the difference between"
					" guidelength and elements of the seglength array is:"
					"%e consider changes the parameters l or seglength \n",
					NAME_CURRENT_COMP,sumOfelements-guideInfo.Length);
	}
	else guideInfo.enableSegments = 0;


///////////////////////////////////////////////////////////////////////////
/////////////// Calculate gravity vector in the guides coordinatesystem
///////////////////////////////////////////////////////////////////////////

/*
	Sets the local gravity vector equal to the global gravity vector (0,-g,0) 
	and when apply the same rotation matrix as applied to guide. 
*/ 
	if (enableGravity != 0){
		Gx0=0, Gy0=-GRAVITY*enableGravity, Gz0=0;
		Coords mcLocG;
		mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,Gy0,0));
		coords_get(mcLocG, &Gx0, &Gy0, &Gz0);
	}
	Circ=2*PI*curvature;
}
#line 12746 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component EndOfelement_5. */
  SIG_MESSAGE("EndOfelement_5 (Init)");

  /* Initializations for component EndOfelement_4. */
  SIG_MESSAGE("EndOfelement_4 (Init)");

  /* Initializations for component elliptical_guide_gravity3. */
  SIG_MESSAGE("elliptical_guide_gravity3 (Init)");
#define mccompcurname  elliptical_guide_gravity3
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 7
#define mvaluesright mccelliptical_guide_gravity3_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity3_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity3_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity3_mvaluesbottom
#define seglength mccelliptical_guide_gravity3_seglength
#define guideInfo mccelliptical_guide_gravity3_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity3_latestParticleCollision
#define Gx mccelliptical_guide_gravity3_Gx
#define Gy mccelliptical_guide_gravity3_Gy
#define Gz mccelliptical_guide_gravity3_Gz
#define Gx0 mccelliptical_guide_gravity3_Gx0
#define Gy0 mccelliptical_guide_gravity3_Gy0
#define Gz0 mccelliptical_guide_gravity3_Gz0
#define Circ mccelliptical_guide_gravity3_Circ
#define dynamicalSegLength mccelliptical_guide_gravity3_dynamicalSegLength
#define l mccelliptical_guide_gravity3_l
#define xwidth mccelliptical_guide_gravity3_xwidth
#define yheight mccelliptical_guide_gravity3_yheight
#define linxw mccelliptical_guide_gravity3_linxw
#define loutxw mccelliptical_guide_gravity3_loutxw
#define linyh mccelliptical_guide_gravity3_linyh
#define loutyh mccelliptical_guide_gravity3_loutyh
#define majorAxisxw mccelliptical_guide_gravity3_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity3_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity3_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity3_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity3_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity3_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity3_dimensionsAt
#define option mccelliptical_guide_gravity3_option
#define R0 mccelliptical_guide_gravity3_R0
#define Qc mccelliptical_guide_gravity3_Qc
#define alpha mccelliptical_guide_gravity3_alpha
#define m mccelliptical_guide_gravity3_m
#define W mccelliptical_guide_gravity3_W
#define alpharight mccelliptical_guide_gravity3_alpharight
#define mright mccelliptical_guide_gravity3_mright
#define alphaleft mccelliptical_guide_gravity3_alphaleft
#define mleft mccelliptical_guide_gravity3_mleft
#define alphatop mccelliptical_guide_gravity3_alphatop
#define mtop mccelliptical_guide_gravity3_mtop
#define alphabottom mccelliptical_guide_gravity3_alphabottom
#define mbottom mccelliptical_guide_gravity3_mbottom
#define verbose mccelliptical_guide_gravity3_verbose
#define enableGravity mccelliptical_guide_gravity3_enableGravity
#define curvature mccelliptical_guide_gravity3_curvature
#line 817 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
	///////////////////////////////////////////////////////////////////////////
	/////////////// Test user input
	///////////////////////////////////////////////////////////////////////////

	if(strcmp(verbose,"on") == 0)
		guideInfo.verboseSetting = 1;
	else guideInfo.verboseSetting = 0;
	

	guideInfo.R0Arr[RightSide] 		= R0;
	guideInfo.QcArr[RightSide] 		= Qc;
	guideInfo.alphaArr[RightSide] 	= alpharight;
	guideInfo.mArr[RightSide] 		= mright;
	guideInfo.WArr[RightSide] 		= W; 

	guideInfo.R0Arr[TopSide] 	 = R0; 
	guideInfo.QcArr[TopSide] 	 = Qc; 
	guideInfo.alphaArr[TopSide]  = alphatop; 
	guideInfo.mArr[TopSide] 	 = mtop;
	guideInfo.WArr[TopSide] 	 = W; 

	guideInfo.R0Arr[LeftSide] 		= R0; 
	guideInfo.QcArr[LeftSide] 		= Qc; 
	guideInfo.alphaArr[LeftSide] 	= alphaleft; 
	guideInfo.mArr[LeftSide] 		= mleft; 
	guideInfo.WArr[LeftSide] 		= W; 

	guideInfo.R0Arr[BottomSide] 	= R0; 
	guideInfo.QcArr[BottomSide] 	= Qc; 
	guideInfo.alphaArr[BottomSide] 	= alphabottom; 
	guideInfo.mArr[BottomSide] 		= mbottom; 
	guideInfo.WArr[BottomSide] 		= W; 

	int sides;
	for (sides = RightSide; sides <= BottomSide; sides++){
		if (guideInfo.alphaArr[sides] == -1) guideInfo.alphaArr[sides] = alpha;
		if (guideInfo.mArr[sides] == -1)	 guideInfo.mArr[sides] = m;
	}
	
	// Test user input for illegal values
	int inputErrors = 0;
	// Lower or equal to zero
	if(l <= 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("length",guideInfo.verboseSetting);
	if(guideInfo.alphaArr[TopSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphatop",guideInfo.verboseSetting);
	if(guideInfo.mArr[TopSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mtop",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[BottomSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphabottom",guideInfo.verboseSetting);
	if(guideInfo.mArr[BottomSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mbottom",guideInfo.verboseSetting);

	if(guideInfo.alphaArr[RightSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alpharight",guideInfo.verboseSetting);
	if(guideInfo.mArr[RightSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mright",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[LeftSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphaleft",guideInfo.verboseSetting);
	if(guideInfo.mArr[LeftSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mleft",guideInfo.verboseSetting);

	// Focal points result in hyperbola instead of an ellipse
	if(l <= -linxw-loutxw) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linw","loutw",linxw,loutxw,guideInfo.verboseSetting);
	if(l <= -linyh-loutyh) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linh","louth",linyh,loutyh,guideInfo.verboseSetting);

	if( strcmp(dimensionsAt,"entrance") != 0 
	 && strcmp(dimensionsAt,"mid") != 0 
	 && strcmp(dimensionsAt,"exit") != 0){
	 	inputErrors += 1;
		printf("dimensionsAt were given an incorrect input."
			   "Input must be string containing \"entrance\",\"mid\" or \"exit\" \n");
	}


	// Terminate program if any input errors occurred
	if(inputErrors != 0 ){
		exit(printf("\nCRITICAL ERROR(S) IN COMPONENT %s" 
				" CONSIDER CHECKING USER INPUT AS %d INPUT ERRORS WAS FOUND.\n",
				NAME_CURRENT_COMP,inputErrors) );
	}


	///////////////////////////////////////////////////////////////////////////
	/////////////// Calculate intern guide values from user input
	///////////////////////////////////////////////////////////////////////////

	/*  Calculate the foci line for the ellipses.
		These can be used to calculate the axes of the ellipses
		using pyth and defination of the ellipse that says distance 
		between the foci and every point on the ellipse is constant.
	*/
	int directDefination = 0;

	if( majorAxisyh != 0 || minorAxisyh != 0 
	 || majorAxisxw != 0 || minorAxisxw != 0)
	{
		directDefination = 1;
		guideInfo.Length = l;	
		
		guideInfo.ellipseMajorAxis[RightSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[RightSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[RightSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[TopSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[TopSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[TopSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[LeftSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[LeftSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[LeftSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[BottomSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[BottomSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[BottomSide] = 0;

		guideInfo.entranceHorizontalWidth =
			2*sqrt(1 - (majorAxisoffsetyh*majorAxisoffsetyh)
			/(majorAxisyh*majorAxisyh) )*minorAxisyh;
		guideInfo.entranceVerticalWidth = 	
			2*sqrt(1 - (majorAxisoffsetxw*majorAxisoffsetxw)
			/(majorAxisxw*majorAxisxw) )*minorAxisxw;
	}

	if ( strcmp(option,"ellipse") == 0 && directDefination == 0)
	{
		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( linxw*linxw + xwidth*xwidth*0.25)
				+ sqrt( (l + loutxw)*(l + loutxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( linyh*linyh + yheight*yheight*0.25)
				+ sqrt( (l + loutyh)*(l + loutyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"exit") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( loutxw*loutxw + xwidth*xwidth*0.25)
				+ sqrt( (l + linxw)*(l + linxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( loutyh*loutyh + yheight*yheight*0.25)
				+ sqrt( (l + linyh)*(l + linyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] =(l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"mid") == 0 ){ 

			guideInfo.Length = l;
			
			guideInfo.ellipseMajorAxis[RightSide] = 
				sqrt( (linxw+l+loutxw)*(linxw+l+loutxw)/4+xwidth*xwidth/4);
			guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = 
				sqrt( (linyh+l+loutyh)*(linyh+l+loutyh)/4+yheight*yheight/4);
			guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];

		}
	}

	guideInfo.entranceHorizontalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[RightSide]
			*guideInfo.ellipseMajorOffset[RightSide]
			/(guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] ) )
		*guideInfo.ellipseMinorAxis[RightSide];
	guideInfo.entranceVerticalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[TopSide]
			*guideInfo.ellipseMajorOffset[TopSide]
			/(guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] ) )
		*guideInfo.ellipseMinorAxis[TopSide];


	if ( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		exit( printf("Critical error in %s; the option for option = halfellipse is currently disabled.",NAME_CURRENT_COMP) );

		double used_focal_vertical;
		double used_focal_horizontal;
		double major_offset_horizontal = 0;
		double major_offset_vertical = 0;

		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+linyh)*(l+linyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+linxw)*(l+linxw) ); 
			major_offset_vertical = l;	
			major_offset_horizontal = l;
		}
		else {	
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+loutyh)*(l+loutyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+loutxw)*(l+loutxw) ); 
		}

		guideInfo.Length = l;

		guideInfo.ellipseMajorAxis[RightSide] = used_focal_horizontal;
		guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

		guideInfo.ellipseMajorOffset[RightSide] = major_offset_horizontal;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = 
								guideInfo.ellipseMajorAxis[RightSide];
		guideInfo.ellipseMinorAxis[LeftSide] = 
								guideInfo.ellipseMinorAxis[RightSide];
		guideInfo.ellipseMajorOffset[LeftSide] = 
								guideInfo.ellipseMajorOffset[RightSide];
		guideInfo.ellipseMinorOffset[LeftSide] = 
								guideInfo.ellipseMinorOffset[RightSide];

		guideInfo.ellipseMajorAxis[TopSide] = used_focal_vertical;
		guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

		guideInfo.ellipseMajorOffset[TopSide] = major_offset_vertical;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = 
								guideInfo.ellipseMajorAxis[TopSide];
		guideInfo.ellipseMinorAxis[BottomSide] = 
								guideInfo.ellipseMinorAxis[TopSide];
		guideInfo.ellipseMajorOffset[BottomSide] = 
								guideInfo.ellipseMajorOffset[TopSide];
		guideInfo.ellipseMinorOffset[BottomSide] = 
								guideInfo.ellipseMinorOffset[TopSide];
	}

	// Applies the properties of the mirrors in the guide given by the user.
	// These variables are used in the reflection functions.

	
	// Sets the mirror type of the guides mirrors
	// These variables are used in the collision functions 
	// to find the type of collision

	// guideInfo.OuterSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[BottomSide] = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_bottom_side_mirror,guideInfo.verboseSetting);

	// guideInfo.InnerSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[BottomSide] =	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_bottom_side_mirror,guideInfo.verboseSetting);

	// Give a warning if all side of the guide is turned off, 
	// as the guide is essentially turned off
	if(    guideInfo.OuterSide[RightSide] 	== 1 
		&& guideInfo.OuterSide[TopSide] 	== 1
		&& guideInfo.OuterSide[LeftSide] 	== 1 
		&& guideInfo.OuterSide[BottomSide]	== 1
		&& guideInfo.InnerSide[RightSide] 	== 1 
		&& guideInfo.InnerSide[TopSide] 	== 1
		&& guideInfo.InnerSide[LeftSide] 	== 1 
		&& guideInfo.InnerSide[BottomSide]  == 1)
		printf("Warning: In %s all the sides of the guide has been disabled,"
				" so it not possible for any particle"
				" to collide with the guide, consider"
				" disabling this component",NAME_CURRENT_COMP);
	
	if(guideInfo.mArr[RightSide] <= 0)  guideInfo.InnerSide[RightSide] = 
														MirrorTypeabsorption;
	if(guideInfo.mArr[TopSide] <= 0) 	 guideInfo.InnerSide[TopSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[LeftSide] <= 0) 	 guideInfo.InnerSide[LeftSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[BottomSide] <= 0) guideInfo.InnerSide[BottomSide] = 
														MirrorTypeabsorption;
	/* if(directDefination == 0){ */
	/* 	guideInfo.entranceHorizontalWidth = xwidth; */
	/* 	guideInfo.entranceVerticalWidth = yheight; */
	/* } */
	
	if( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		guideInfo.entranceHorizontalWidth =
			(guideInfo.ellipseMinorAxis[RightSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[RightSide]
						*guideInfo.ellipseMajorOffset[RightSide] )
			/( guideInfo.ellipseMajorAxis[RightSide]
			  *guideInfo.ellipseMajorAxis[RightSide] ) )
			+ guideInfo.ellipseMinorOffset[RightSide] )*2;
		guideInfo.entranceVerticalWidth = 
			(guideInfo.ellipseMinorAxis[TopSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[TopSide]
						*guideInfo.ellipseMajorOffset[TopSide] )
			/( guideInfo.ellipseMajorAxis[TopSide]
			  *guideInfo.ellipseMajorAxis[TopSide] ) )
			+ guideInfo.ellipseMinorOffset[TopSide] )*2;
	}


	guideInfo.EnclosingBoxOn = 0;

	/*
	double DefaultArray1[8] = { 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0};
	double DefaultArray2[8] = { 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0};
	double DefaultArray3[8] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0}; 

	guideInfo.EnclosingBoxOn = 0;
	double *xinput;
	if ( xInput != NULL ){ xinput = xInput; guideInfo.EnclosingBoxOn = 1; }
	else { xinput = DefaultArray1; }
	double *yinput;
	if ( yInput != NULL ){ yinput = yInput; guideInfo.EnclosingBoxOn = 1;}
	else { yinput = DefaultArray2; }
	double *zinput;
	if ( zInput != NULL ){ zinput = zInput; guideInfo.EnclosingBoxOn = 1;}
	else { zinput = DefaultArray3; }
	*/

	/*
	double xarray[8] ={ guideInfo.ellipseMinorAxis[0]*xinput[0],
						guideInfo.ellipseMinorAxis[2]*xinput[1],
						guideInfo.ellipseMinorAxis[2]*xinput[2],
						guideInfo.ellipseMinorAxis[0]*xinput[3],
						guideInfo.ellipseMinorAxis[0]*xinput[4],
						guideInfo.ellipseMinorAxis[2]*xinput[5],
						guideInfo.ellipseMinorAxis[2]*xinput[6],
						guideInfo.ellipseMinorAxis[0]*xinput[7] };
	double yarray[8] ={ guideInfo.ellipseMinorAxis[1]*yinput[0],
						guideInfo.ellipseMinorAxis[1]*yinput[1],
						guideInfo.ellipseMinorAxis[3]*yinput[2],
						guideInfo.ellipseMinorAxis[3]*yinput[3],
						guideInfo.ellipseMinorAxis[1]*yinput[4],
						guideInfo.ellipseMinorAxis[1]*yinput[5],
						guideInfo.ellipseMinorAxis[3]*yinput[6],
						guideInfo.ellipseMinorAxis[3]*yinput[7] };
	double zarray[8] ={ guideInfo.Length/2*zinput[0]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[1]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[2]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[3]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[4]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[5]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[6]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[7]+guideInfo.Length/2 };
	int i = 0;
	for(i = 0; i < 8; i++){
		guideInfo.xArray[i] = xarray[i];
		guideInfo.yArray[i] = yarray[i];
		guideInfo.zArray[i] = zarray[i];
	}
	*/

	guideInfo.exitVerticalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide]) 
					)/(guideInfo.ellipseMajorAxis[BottomSide] 
					*guideInfo.ellipseMajorAxis[BottomSide] )
		)*guideInfo.ellipseMinorAxis[BottomSide];

	guideInfo.exitHorizontalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide]) 
					)/(guideInfo.ellipseMajorAxis[RightSide] 
					*guideInfo.ellipseMajorAxis[RightSide] ) 
		)*guideInfo.ellipseMinorAxis[RightSide];

    //////////////////segmentation of m values

	// Are the arrays empty?
	if(mvaluesright != NULL || mvaluesleft != NULL 
	|| mvaluestop != NULL || mvaluesbottom != NULL)
	{
		guideInfo.enableSegments = 1;

		guideInfo.numberOfSegments = sizeof(mvaluesright)/sizeof(mvaluesright[0]);

		//printf("Length is %i\n",guideInfo.numberOfSegments);
		guideInfo.mValuesright = mvaluesright;
		guideInfo.mValuesleft = mvaluesleft;
		guideInfo.mValuestop = mvaluestop;
		guideInfo.mValuesbottom = mvaluesbottom;
		//printf("Seglength ... %f %f %f\n",seglength[0],seglength[1],seglength[2]);

		// Are the arrays of equal length?
		if(seglength == NULL){
			dynamicalSegLength = 
					realloc(dynamicalSegLength,
							guideInfo.numberOfSegments*sizeof(double) );
			int i;
			for (i = 0; i < guideInfo.numberOfSegments; ++i){
				dynamicalSegLength[i] = 
						guideInfo.Length/guideInfo.numberOfSegments;
				}
			guideInfo.segLength = dynamicalSegLength;
		}
		else guideInfo.segLength = seglength;

		if( guideInfo.numberOfSegments != sizeof(mvaluesright)/sizeof(mvaluesright[0]) 
		 || guideInfo.numberOfSegments != sizeof(mvaluesleft)/sizeof(mvaluesleft[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluestop)/sizeof(mvaluestop[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0])
		 || (guideInfo.segLength == NULL
		  & guideInfo.numberOfSegments != sizeof(seglength)/sizeof(guideInfo.segLength[0]) 
		  ) ) {
			
			printf("Error in userinput inside %s, the length of the arrays"
				   " mvalues and seglength are not equal\n",NAME_CURRENT_COMP);
			printf("The length of the arrays are: mValuesright is %lu,"
				   " mvaluesleft is %lu, mvaluestop is %lu, mvaluesbottom is"
				   " %lu and seglength is %lu and should be %d \n; Above assume that the arrays are using double \n",
					sizeof(mvaluesright)/sizeof(double),
					sizeof(mvaluesleft)/sizeof(double),
					sizeof(mvaluestop)/sizeof(double),
					sizeof(mvaluesbottom)/sizeof(double),
					sizeof(guideInfo.segLength)/sizeof(double),
					guideInfo.numberOfSegments	
					);
				
			if ( guideInfo.verboseSetting )	{
				int i;
				
				printf("The Values of mvaluesright is: [");
				for(i=0; i < sizeof(mvaluesright)/sizeof(mvaluesright[0]); i++) 	
						printf("%e,",guideInfo.mValuesright[i] );
				printf("]\n");
				
				printf("The Values of mvaluesleft is: [");
				for(i=0; i < sizeof(mvaluesleft)/sizeof(mvaluesleft[0]); i++) 	
						printf("%e,",guideInfo.mValuesleft[i] );
				printf("]\n");
				
				printf("The Values of mvaluestop is: [");
				for(i=0; i < sizeof(mvaluestop)/sizeof(mvaluestop[0]); i++) 	
						printf("%e,",guideInfo.mValuestop[i] );
				printf("]\n");

				printf("The Values of mvaluesbottom is: [");
				for(i=0; i < sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0]); i++) 
						printf("%e,",guideInfo.mValuesbottom[i] );
				printf("]\n");

				printf("The Values of seglength is: [");
				for(i=0; i < sizeof(guideInfo.segLength)/sizeof(guideInfo.segLength[0]); i++) 
						printf("%e,",guideInfo.segLength[i]);
				printf("]\n");
			}
			exit( printf("Exit due to critical error in userinput for the"
						" component %s, consider having a look at the input"
						" for following: mvaluesright,mvaluesleft,mvaluestop,"
						"mvaluesbottom and/or seglength.",NAME_CURRENT_COMP) );
		}
		//
		double sumOfelements=0;
		int i;
		for(i=0;i< guideInfo.numberOfSegments; i++) {
			sumOfelements += guideInfo.segLength[i];
		}
		if ( guideInfo.verboseSetting 
		&& fabs(sumOfelements-guideInfo.Length) > 1e-9 )
			printf("Error in userinput inside %s, the difference between"
					" guidelength and elements of the seglength array is:"
					"%e consider changes the parameters l or seglength \n",
					NAME_CURRENT_COMP,sumOfelements-guideInfo.Length);
	}
	else guideInfo.enableSegments = 0;


///////////////////////////////////////////////////////////////////////////
/////////////// Calculate gravity vector in the guides coordinatesystem
///////////////////////////////////////////////////////////////////////////

/*
	Sets the local gravity vector equal to the global gravity vector (0,-g,0) 
	and when apply the same rotation matrix as applied to guide. 
*/ 
	if (enableGravity != 0){
		Gx0=0, Gy0=-GRAVITY*enableGravity, Gz0=0;
		Coords mcLocG;
		mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,Gy0,0));
		coords_get(mcLocG, &Gx0, &Gy0, &Gz0);
	}
	Circ=2*PI*curvature;
}
#line 13446 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component EndOfelement_3. */
  SIG_MESSAGE("EndOfelement_3 (Init)");

  /* Initializations for component Curved_guide_2. */
  SIG_MESSAGE("Curved_guide_2 (Init)");
#define mccompcurname  Curved_guide_2
#define mccompcurtype  Bender
#define mccompcurindex 9
#define bk mccCurved_guide_2_bk
#define mWin mccCurved_guide_2_mWin
#define w mccCurved_guide_2_w
#define h mccCurved_guide_2_h
#define r mccCurved_guide_2_r
#define Win mccCurved_guide_2_Win
#define k mccCurved_guide_2_k
#define d mccCurved_guide_2_d
#define l mccCurved_guide_2_l
#define R0a mccCurved_guide_2_R0a
#define Qca mccCurved_guide_2_Qca
#define alphaa mccCurved_guide_2_alphaa
#define ma mccCurved_guide_2_ma
#define Wa mccCurved_guide_2_Wa
#define R0i mccCurved_guide_2_R0i
#define Qci mccCurved_guide_2_Qci
#define alphai mccCurved_guide_2_alphai
#define mi mccCurved_guide_2_mi
#define Wi mccCurved_guide_2_Wi
#define R0s mccCurved_guide_2_R0s
#define Qcs mccCurved_guide_2_Qcs
#define alphas mccCurved_guide_2_alphas
#define ms mccCurved_guide_2_ms
#define Ws mccCurved_guide_2_Ws
#line 112 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
if (r <0)
      { fprintf(stderr,"Bender: error: %s: to bend in the other direction\n", NAME_CURRENT_COMP);
        fprintf(stderr,"        rotate comp on z-axis by 180 deg.\n"); exit(-1); }

      if (k*d > w)
      { fprintf(stderr,"Bender: error: %s has (k*d > w).\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (w*h*r*Win*k == 0)
      { fprintf(stderr,"Bender: error: %s has one of w,h,r,Win,k null.\n", NAME_CURRENT_COMP);
        exit(-1); }
      /* width of one channel + thickness d of partition */
      mWin = Win;
      if (l!= 0 && r != 0) mWin = (double)l/(double)r;
      bk=(w+d)/k;
      if (mcgravitation) fprintf(stderr,"WARNING: Bender: %s: "
        "This component produces wrong results with gravitation !\n",
        NAME_CURRENT_COMP);
}
#line 13549 "./FullInstrument_v4.7_frontend.c"
#undef Ws
#undef ms
#undef alphas
#undef Qcs
#undef R0s
#undef Wi
#undef mi
#undef alphai
#undef Qci
#undef R0i
#undef Wa
#undef ma
#undef alphaa
#undef Qca
#undef R0a
#undef l
#undef d
#undef k
#undef Win
#undef r
#undef h
#undef w
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component EndOfelement_2. */
  SIG_MESSAGE("EndOfelement_2 (Init)");

  /* Initializations for component elliptical_guide_gravity1. */
  SIG_MESSAGE("elliptical_guide_gravity1 (Init)");
#define mccompcurname  elliptical_guide_gravity1
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 11
#define mvaluesright mccelliptical_guide_gravity1_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity1_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity1_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity1_mvaluesbottom
#define seglength mccelliptical_guide_gravity1_seglength
#define guideInfo mccelliptical_guide_gravity1_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity1_latestParticleCollision
#define Gx mccelliptical_guide_gravity1_Gx
#define Gy mccelliptical_guide_gravity1_Gy
#define Gz mccelliptical_guide_gravity1_Gz
#define Gx0 mccelliptical_guide_gravity1_Gx0
#define Gy0 mccelliptical_guide_gravity1_Gy0
#define Gz0 mccelliptical_guide_gravity1_Gz0
#define Circ mccelliptical_guide_gravity1_Circ
#define dynamicalSegLength mccelliptical_guide_gravity1_dynamicalSegLength
#define l mccelliptical_guide_gravity1_l
#define xwidth mccelliptical_guide_gravity1_xwidth
#define yheight mccelliptical_guide_gravity1_yheight
#define linxw mccelliptical_guide_gravity1_linxw
#define loutxw mccelliptical_guide_gravity1_loutxw
#define linyh mccelliptical_guide_gravity1_linyh
#define loutyh mccelliptical_guide_gravity1_loutyh
#define majorAxisxw mccelliptical_guide_gravity1_majorAxisxw
#define minorAxisxw mccelliptical_guide_gravity1_minorAxisxw
#define majorAxisyh mccelliptical_guide_gravity1_majorAxisyh
#define minorAxisyh mccelliptical_guide_gravity1_minorAxisyh
#define majorAxisoffsetxw mccelliptical_guide_gravity1_majorAxisoffsetxw
#define majorAxisoffsetyh mccelliptical_guide_gravity1_majorAxisoffsetyh
#define dimensionsAt mccelliptical_guide_gravity1_dimensionsAt
#define option mccelliptical_guide_gravity1_option
#define R0 mccelliptical_guide_gravity1_R0
#define Qc mccelliptical_guide_gravity1_Qc
#define alpha mccelliptical_guide_gravity1_alpha
#define m mccelliptical_guide_gravity1_m
#define W mccelliptical_guide_gravity1_W
#define alpharight mccelliptical_guide_gravity1_alpharight
#define mright mccelliptical_guide_gravity1_mright
#define alphaleft mccelliptical_guide_gravity1_alphaleft
#define mleft mccelliptical_guide_gravity1_mleft
#define alphatop mccelliptical_guide_gravity1_alphatop
#define mtop mccelliptical_guide_gravity1_mtop
#define alphabottom mccelliptical_guide_gravity1_alphabottom
#define mbottom mccelliptical_guide_gravity1_mbottom
#define verbose mccelliptical_guide_gravity1_verbose
#define enableGravity mccelliptical_guide_gravity1_enableGravity
#define curvature mccelliptical_guide_gravity1_curvature
#line 817 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
	///////////////////////////////////////////////////////////////////////////
	/////////////// Test user input
	///////////////////////////////////////////////////////////////////////////

	if(strcmp(verbose,"on") == 0)
		guideInfo.verboseSetting = 1;
	else guideInfo.verboseSetting = 0;
	

	guideInfo.R0Arr[RightSide] 		= R0;
	guideInfo.QcArr[RightSide] 		= Qc;
	guideInfo.alphaArr[RightSide] 	= alpharight;
	guideInfo.mArr[RightSide] 		= mright;
	guideInfo.WArr[RightSide] 		= W; 

	guideInfo.R0Arr[TopSide] 	 = R0; 
	guideInfo.QcArr[TopSide] 	 = Qc; 
	guideInfo.alphaArr[TopSide]  = alphatop; 
	guideInfo.mArr[TopSide] 	 = mtop;
	guideInfo.WArr[TopSide] 	 = W; 

	guideInfo.R0Arr[LeftSide] 		= R0; 
	guideInfo.QcArr[LeftSide] 		= Qc; 
	guideInfo.alphaArr[LeftSide] 	= alphaleft; 
	guideInfo.mArr[LeftSide] 		= mleft; 
	guideInfo.WArr[LeftSide] 		= W; 

	guideInfo.R0Arr[BottomSide] 	= R0; 
	guideInfo.QcArr[BottomSide] 	= Qc; 
	guideInfo.alphaArr[BottomSide] 	= alphabottom; 
	guideInfo.mArr[BottomSide] 		= mbottom; 
	guideInfo.WArr[BottomSide] 		= W; 

	int sides;
	for (sides = RightSide; sides <= BottomSide; sides++){
		if (guideInfo.alphaArr[sides] == -1) guideInfo.alphaArr[sides] = alpha;
		if (guideInfo.mArr[sides] == -1)	 guideInfo.mArr[sides] = m;
	}
	
	// Test user input for illegal values
	int inputErrors = 0;
	// Lower or equal to zero
	if(l <= 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("length",guideInfo.verboseSetting);
	if(guideInfo.alphaArr[TopSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphatop",guideInfo.verboseSetting);
	if(guideInfo.mArr[TopSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mtop",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[BottomSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphabottom",guideInfo.verboseSetting);
	if(guideInfo.mArr[BottomSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mbottom",guideInfo.verboseSetting);

	if(guideInfo.alphaArr[RightSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alpharight",guideInfo.verboseSetting);
	if(guideInfo.mArr[RightSide] < 0)       inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mright",guideInfo.verboseSetting);
	
	if(guideInfo.alphaArr[LeftSide] < 0) inputErrors += 
				guide_elliptical_illegalInputLessThanZero("alphaleft",guideInfo.verboseSetting);
	if(guideInfo.mArr[LeftSide] < 0)      inputErrors += 
				guide_elliptical_illegalInputLessThanZero("mleft",guideInfo.verboseSetting);

	// Focal points result in hyperbola instead of an ellipse
	if(l <= -linxw-loutxw) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linw","loutw",linxw,loutxw,guideInfo.verboseSetting);
	if(l <= -linyh-loutyh) inputErrors += guide_elliptical_illegalInputFocalPointsHyperbola(
						"linh","louth",linyh,loutyh,guideInfo.verboseSetting);

	if( strcmp(dimensionsAt,"entrance") != 0 
	 && strcmp(dimensionsAt,"mid") != 0 
	 && strcmp(dimensionsAt,"exit") != 0){
	 	inputErrors += 1;
		printf("dimensionsAt were given an incorrect input."
			   "Input must be string containing \"entrance\",\"mid\" or \"exit\" \n");
	}


	// Terminate program if any input errors occurred
	if(inputErrors != 0 ){
		exit(printf("\nCRITICAL ERROR(S) IN COMPONENT %s" 
				" CONSIDER CHECKING USER INPUT AS %d INPUT ERRORS WAS FOUND.\n",
				NAME_CURRENT_COMP,inputErrors) );
	}


	///////////////////////////////////////////////////////////////////////////
	/////////////// Calculate intern guide values from user input
	///////////////////////////////////////////////////////////////////////////

	/*  Calculate the foci line for the ellipses.
		These can be used to calculate the axes of the ellipses
		using pyth and defination of the ellipse that says distance 
		between the foci and every point on the ellipse is constant.
	*/
	int directDefination = 0;

	if( majorAxisyh != 0 || minorAxisyh != 0 
	 || majorAxisxw != 0 || minorAxisxw != 0)
	{
		directDefination = 1;
		guideInfo.Length = l;	
		
		guideInfo.ellipseMajorAxis[RightSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[RightSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[RightSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[TopSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[TopSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[TopSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = majorAxisxw;
		guideInfo.ellipseMinorAxis[LeftSide] = minorAxisxw;
		guideInfo.ellipseMajorOffset[LeftSide] = majorAxisoffsetxw;
		guideInfo.ellipseMinorOffset[LeftSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = majorAxisyh;
		guideInfo.ellipseMinorAxis[BottomSide] = minorAxisyh;
		guideInfo.ellipseMajorOffset[BottomSide] = majorAxisoffsetyh;
		guideInfo.ellipseMinorOffset[BottomSide] = 0;

		guideInfo.entranceHorizontalWidth =
			2*sqrt(1 - (majorAxisoffsetyh*majorAxisoffsetyh)
			/(majorAxisyh*majorAxisyh) )*minorAxisyh;
		guideInfo.entranceVerticalWidth = 	
			2*sqrt(1 - (majorAxisoffsetxw*majorAxisoffsetxw)
			/(majorAxisxw*majorAxisxw) )*minorAxisxw;
	}

	if ( strcmp(option,"ellipse") == 0 && directDefination == 0)
	{
		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( linxw*linxw + xwidth*xwidth*0.25)
				+ sqrt( (l + loutxw)*(l + loutxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( linyh*linyh + yheight*yheight*0.25)
				+ sqrt( (l + loutyh)*(l + loutyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"exit") == 0 ){ 
			double lofbs_horizontal = 
				  sqrt( loutxw*loutxw + xwidth*xwidth*0.25)
				+ sqrt( (l + linxw)*(l + linxw) + xwidth*xwidth*0.25);

			double lofbs_vertical = 
				  sqrt( loutyh*loutyh + yheight*yheight*0.25)
				+ sqrt( (l + linyh)*(l + linyh) + yheight*yheight*0.25);

			guideInfo.Length = l;

			guideInfo.ellipseMajorAxis[RightSide] = lofbs_horizontal/2;
			guideInfo.ellipseMinorAxis[RightSide] =
				sqrt(0.25*lofbs_horizontal*lofbs_horizontal
					-0.25*(l+linxw+loutxw)*(l+linxw+loutxw) );

			guideInfo.ellipseMajorOffset[RightSide] =(l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = lofbs_vertical/2;
			
			guideInfo.ellipseMinorAxis[TopSide] = 
					sqrt(0.25*lofbs_vertical*lofbs_vertical
						-0.25*(l+linyh+loutyh)*(l+linyh+loutyh)	);

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];
		}
		if ( strcmp(dimensionsAt,"mid") == 0 ){ 

			guideInfo.Length = l;
			
			guideInfo.ellipseMajorAxis[RightSide] = 
				sqrt( (linxw+l+loutxw)*(linxw+l+loutxw)/4+xwidth*xwidth/4);
			guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

			guideInfo.ellipseMajorOffset[RightSide] = (l+linxw+loutxw)/2-linxw;
			guideInfo.ellipseMinorOffset[RightSide] = 0;

			guideInfo.ellipseMajorAxis[LeftSide] = 
										guideInfo.ellipseMajorAxis[RightSide];
			guideInfo.ellipseMinorAxis[LeftSide] = 
										guideInfo.ellipseMinorAxis[RightSide];
			guideInfo.ellipseMajorOffset[LeftSide] = 
										guideInfo.ellipseMajorOffset[RightSide];
			guideInfo.ellipseMinorOffset[LeftSide] = 
										guideInfo.ellipseMinorOffset[RightSide];

			guideInfo.ellipseMajorAxis[TopSide] = 
				sqrt( (linyh+l+loutyh)*(linyh+l+loutyh)/4+yheight*yheight/4);
			guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

			guideInfo.ellipseMajorOffset[TopSide] = (l+linyh+loutyh)/2-linyh;
			guideInfo.ellipseMinorOffset[TopSide] = 0;

			guideInfo.ellipseMajorAxis[BottomSide] = 
										guideInfo.ellipseMajorAxis[TopSide];
			guideInfo.ellipseMinorAxis[BottomSide] = 
										guideInfo.ellipseMinorAxis[TopSide];
			guideInfo.ellipseMajorOffset[BottomSide] = 
										guideInfo.ellipseMajorOffset[TopSide];
			guideInfo.ellipseMinorOffset[BottomSide] = 
										guideInfo.ellipseMinorOffset[TopSide];

		}
	}

	guideInfo.entranceHorizontalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[RightSide]
			*guideInfo.ellipseMajorOffset[RightSide]
			/(guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] ) )
		*guideInfo.ellipseMinorAxis[RightSide];
	guideInfo.entranceVerticalWidth = 2*sqrt(
		1 - guideInfo.ellipseMajorOffset[TopSide]
			*guideInfo.ellipseMajorOffset[TopSide]
			/(guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] ) )
		*guideInfo.ellipseMinorAxis[TopSide];


	if ( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		exit( printf("Critical error in %s; the option for option = halfellipse is currently disabled.",NAME_CURRENT_COMP) );

		double used_focal_vertical;
		double used_focal_horizontal;
		double major_offset_horizontal = 0;
		double major_offset_vertical = 0;

		if ( strcmp(dimensionsAt,"entrance") == 0 ){ 
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+linyh)*(l+linyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+linxw)*(l+linxw) ); 
			major_offset_vertical = l;	
			major_offset_horizontal = l;
		}
		else {	
			used_focal_vertical = sqrt( (yheight*yheight)/4 
								+ (l+loutyh)*(l+loutyh) );
			used_focal_horizontal = sqrt( (xwidth*xwidth)/4 
								+ (l+loutxw)*(l+loutxw) ); 
		}

		guideInfo.Length = l;

		guideInfo.ellipseMajorAxis[RightSide] = used_focal_horizontal;
		guideInfo.ellipseMinorAxis[RightSide] = xwidth/2;

		guideInfo.ellipseMajorOffset[RightSide] = major_offset_horizontal;
		guideInfo.ellipseMinorOffset[RightSide] = 0;

		guideInfo.ellipseMajorAxis[LeftSide] = 
								guideInfo.ellipseMajorAxis[RightSide];
		guideInfo.ellipseMinorAxis[LeftSide] = 
								guideInfo.ellipseMinorAxis[RightSide];
		guideInfo.ellipseMajorOffset[LeftSide] = 
								guideInfo.ellipseMajorOffset[RightSide];
		guideInfo.ellipseMinorOffset[LeftSide] = 
								guideInfo.ellipseMinorOffset[RightSide];

		guideInfo.ellipseMajorAxis[TopSide] = used_focal_vertical;
		guideInfo.ellipseMinorAxis[TopSide] = yheight/2;

		guideInfo.ellipseMajorOffset[TopSide] = major_offset_vertical;
		guideInfo.ellipseMinorOffset[TopSide] = 0;

		guideInfo.ellipseMajorAxis[BottomSide] = 
								guideInfo.ellipseMajorAxis[TopSide];
		guideInfo.ellipseMinorAxis[BottomSide] = 
								guideInfo.ellipseMinorAxis[TopSide];
		guideInfo.ellipseMajorOffset[BottomSide] = 
								guideInfo.ellipseMajorOffset[TopSide];
		guideInfo.ellipseMinorOffset[BottomSide] = 
								guideInfo.ellipseMinorOffset[TopSide];
	}

	// Applies the properties of the mirrors in the guide given by the user.
	// These variables are used in the reflection functions.

	
	// Sets the mirror type of the guides mirrors
	// These variables are used in the collision functions 
	// to find the type of collision

	// guideInfo.OuterSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.OuterSide[BottomSide] = 	
	// 	guide_elliptical_getMirrorTypeFromInput(outer_bottom_side_mirror,guideInfo.verboseSetting);

	// guideInfo.InnerSide[RightSide]  = 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_right_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[TopSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_top_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[LeftSide]  	= 	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_left_side_mirror,guideInfo.verboseSetting);
	// guideInfo.InnerSide[BottomSide] =	
	// 	guide_elliptical_getMirrorTypeFromInput(inner_bottom_side_mirror,guideInfo.verboseSetting);

	// Give a warning if all side of the guide is turned off, 
	// as the guide is essentially turned off
	if(    guideInfo.OuterSide[RightSide] 	== 1 
		&& guideInfo.OuterSide[TopSide] 	== 1
		&& guideInfo.OuterSide[LeftSide] 	== 1 
		&& guideInfo.OuterSide[BottomSide]	== 1
		&& guideInfo.InnerSide[RightSide] 	== 1 
		&& guideInfo.InnerSide[TopSide] 	== 1
		&& guideInfo.InnerSide[LeftSide] 	== 1 
		&& guideInfo.InnerSide[BottomSide]  == 1)
		printf("Warning: In %s all the sides of the guide has been disabled,"
				" so it not possible for any particle"
				" to collide with the guide, consider"
				" disabling this component",NAME_CURRENT_COMP);
	
	if(guideInfo.mArr[RightSide] <= 0)  guideInfo.InnerSide[RightSide] = 
														MirrorTypeabsorption;
	if(guideInfo.mArr[TopSide] <= 0) 	 guideInfo.InnerSide[TopSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[LeftSide] <= 0) 	 guideInfo.InnerSide[LeftSide] 	= 
														MirrorTypeabsorption;
	if(guideInfo.mArr[BottomSide] <= 0) guideInfo.InnerSide[BottomSide] = 
														MirrorTypeabsorption;
	/* if(directDefination == 0){ */
	/* 	guideInfo.entranceHorizontalWidth = xwidth; */
	/* 	guideInfo.entranceVerticalWidth = yheight; */
	/* } */
	
	if( strcmp(option,"halfellipse") == 0 && directDefination == 0 ){
		guideInfo.entranceHorizontalWidth =
			(guideInfo.ellipseMinorAxis[RightSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[RightSide]
						*guideInfo.ellipseMajorOffset[RightSide] )
			/( guideInfo.ellipseMajorAxis[RightSide]
			  *guideInfo.ellipseMajorAxis[RightSide] ) )
			+ guideInfo.ellipseMinorOffset[RightSide] )*2;
		guideInfo.entranceVerticalWidth = 
			(guideInfo.ellipseMinorAxis[TopSide]
			* sqrt(1 - ( guideInfo.ellipseMajorOffset[TopSide]
						*guideInfo.ellipseMajorOffset[TopSide] )
			/( guideInfo.ellipseMajorAxis[TopSide]
			  *guideInfo.ellipseMajorAxis[TopSide] ) )
			+ guideInfo.ellipseMinorOffset[TopSide] )*2;
	}


	guideInfo.EnclosingBoxOn = 0;

	/*
	double DefaultArray1[8] = { 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0};
	double DefaultArray2[8] = { 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0};
	double DefaultArray3[8] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0}; 

	guideInfo.EnclosingBoxOn = 0;
	double *xinput;
	if ( xInput != NULL ){ xinput = xInput; guideInfo.EnclosingBoxOn = 1; }
	else { xinput = DefaultArray1; }
	double *yinput;
	if ( yInput != NULL ){ yinput = yInput; guideInfo.EnclosingBoxOn = 1;}
	else { yinput = DefaultArray2; }
	double *zinput;
	if ( zInput != NULL ){ zinput = zInput; guideInfo.EnclosingBoxOn = 1;}
	else { zinput = DefaultArray3; }
	*/

	/*
	double xarray[8] ={ guideInfo.ellipseMinorAxis[0]*xinput[0],
						guideInfo.ellipseMinorAxis[2]*xinput[1],
						guideInfo.ellipseMinorAxis[2]*xinput[2],
						guideInfo.ellipseMinorAxis[0]*xinput[3],
						guideInfo.ellipseMinorAxis[0]*xinput[4],
						guideInfo.ellipseMinorAxis[2]*xinput[5],
						guideInfo.ellipseMinorAxis[2]*xinput[6],
						guideInfo.ellipseMinorAxis[0]*xinput[7] };
	double yarray[8] ={ guideInfo.ellipseMinorAxis[1]*yinput[0],
						guideInfo.ellipseMinorAxis[1]*yinput[1],
						guideInfo.ellipseMinorAxis[3]*yinput[2],
						guideInfo.ellipseMinorAxis[3]*yinput[3],
						guideInfo.ellipseMinorAxis[1]*yinput[4],
						guideInfo.ellipseMinorAxis[1]*yinput[5],
						guideInfo.ellipseMinorAxis[3]*yinput[6],
						guideInfo.ellipseMinorAxis[3]*yinput[7] };
	double zarray[8] ={ guideInfo.Length/2*zinput[0]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[1]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[2]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[3]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[4]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[5]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[6]+guideInfo.Length/2,
						guideInfo.Length/2*zinput[7]+guideInfo.Length/2 };
	int i = 0;
	for(i = 0; i < 8; i++){
		guideInfo.xArray[i] = xarray[i];
		guideInfo.yArray[i] = yarray[i];
		guideInfo.zArray[i] = zarray[i];
	}
	*/

	guideInfo.exitVerticalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[BottomSide]) 
					)/(guideInfo.ellipseMajorAxis[BottomSide] 
					*guideInfo.ellipseMajorAxis[BottomSide] )
		)*guideInfo.ellipseMinorAxis[BottomSide];

	guideInfo.exitHorizontalWidth =
		2*sqrt(1 - ( (guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide])
					*(guideInfo.Length-guideInfo.ellipseMajorOffset[RightSide]) 
					)/(guideInfo.ellipseMajorAxis[RightSide] 
					*guideInfo.ellipseMajorAxis[RightSide] ) 
		)*guideInfo.ellipseMinorAxis[RightSide];

    //////////////////segmentation of m values

	// Are the arrays empty?
	if(mvaluesright != NULL || mvaluesleft != NULL 
	|| mvaluestop != NULL || mvaluesbottom != NULL)
	{
		guideInfo.enableSegments = 1;

		guideInfo.numberOfSegments = sizeof(mvaluesright)/sizeof(mvaluesright[0]);

		//printf("Length is %i\n",guideInfo.numberOfSegments);
		guideInfo.mValuesright = mvaluesright;
		guideInfo.mValuesleft = mvaluesleft;
		guideInfo.mValuestop = mvaluestop;
		guideInfo.mValuesbottom = mvaluesbottom;
		//printf("Seglength ... %f %f %f\n",seglength[0],seglength[1],seglength[2]);

		// Are the arrays of equal length?
		if(seglength == NULL){
			dynamicalSegLength = 
					realloc(dynamicalSegLength,
							guideInfo.numberOfSegments*sizeof(double) );
			int i;
			for (i = 0; i < guideInfo.numberOfSegments; ++i){
				dynamicalSegLength[i] = 
						guideInfo.Length/guideInfo.numberOfSegments;
				}
			guideInfo.segLength = dynamicalSegLength;
		}
		else guideInfo.segLength = seglength;

		if( guideInfo.numberOfSegments != sizeof(mvaluesright)/sizeof(mvaluesright[0]) 
		 || guideInfo.numberOfSegments != sizeof(mvaluesleft)/sizeof(mvaluesleft[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluestop)/sizeof(mvaluestop[0])
		 || guideInfo.numberOfSegments != sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0])
		 || (guideInfo.segLength == NULL
		  & guideInfo.numberOfSegments != sizeof(seglength)/sizeof(guideInfo.segLength[0]) 
		  ) ) {
			
			printf("Error in userinput inside %s, the length of the arrays"
				   " mvalues and seglength are not equal\n",NAME_CURRENT_COMP);
			printf("The length of the arrays are: mValuesright is %lu,"
				   " mvaluesleft is %lu, mvaluestop is %lu, mvaluesbottom is"
				   " %lu and seglength is %lu and should be %d \n; Above assume that the arrays are using double \n",
					sizeof(mvaluesright)/sizeof(double),
					sizeof(mvaluesleft)/sizeof(double),
					sizeof(mvaluestop)/sizeof(double),
					sizeof(mvaluesbottom)/sizeof(double),
					sizeof(guideInfo.segLength)/sizeof(double),
					guideInfo.numberOfSegments	
					);
				
			if ( guideInfo.verboseSetting )	{
				int i;
				
				printf("The Values of mvaluesright is: [");
				for(i=0; i < sizeof(mvaluesright)/sizeof(mvaluesright[0]); i++) 	
						printf("%e,",guideInfo.mValuesright[i] );
				printf("]\n");
				
				printf("The Values of mvaluesleft is: [");
				for(i=0; i < sizeof(mvaluesleft)/sizeof(mvaluesleft[0]); i++) 	
						printf("%e,",guideInfo.mValuesleft[i] );
				printf("]\n");
				
				printf("The Values of mvaluestop is: [");
				for(i=0; i < sizeof(mvaluestop)/sizeof(mvaluestop[0]); i++) 	
						printf("%e,",guideInfo.mValuestop[i] );
				printf("]\n");

				printf("The Values of mvaluesbottom is: [");
				for(i=0; i < sizeof(mvaluesbottom)/sizeof(mvaluesbottom[0]); i++) 
						printf("%e,",guideInfo.mValuesbottom[i] );
				printf("]\n");

				printf("The Values of seglength is: [");
				for(i=0; i < sizeof(guideInfo.segLength)/sizeof(guideInfo.segLength[0]); i++) 
						printf("%e,",guideInfo.segLength[i]);
				printf("]\n");
			}
			exit( printf("Exit due to critical error in userinput for the"
						" component %s, consider having a look at the input"
						" for following: mvaluesright,mvaluesleft,mvaluestop,"
						"mvaluesbottom and/or seglength.",NAME_CURRENT_COMP) );
		}
		//
		double sumOfelements=0;
		int i;
		for(i=0;i< guideInfo.numberOfSegments; i++) {
			sumOfelements += guideInfo.segLength[i];
		}
		if ( guideInfo.verboseSetting 
		&& fabs(sumOfelements-guideInfo.Length) > 1e-9 )
			printf("Error in userinput inside %s, the difference between"
					" guidelength and elements of the seglength array is:"
					"%e consider changes the parameters l or seglength \n",
					NAME_CURRENT_COMP,sumOfelements-guideInfo.Length);
	}
	else guideInfo.enableSegments = 0;


///////////////////////////////////////////////////////////////////////////
/////////////// Calculate gravity vector in the guides coordinatesystem
///////////////////////////////////////////////////////////////////////////

/*
	Sets the local gravity vector equal to the global gravity vector (0,-g,0) 
	and when apply the same rotation matrix as applied to guide. 
*/ 
	if (enableGravity != 0){
		Gx0=0, Gy0=-GRAVITY*enableGravity, Gz0=0;
		Coords mcLocG;
		mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,Gy0,0));
		coords_get(mcLocG, &Gx0, &Gy0, &Gz0);
	}
	Circ=2*PI*curvature;
}
#line 14224 "./FullInstrument_v4.7_frontend.c"
#undef curvature
#undef enableGravity
#undef verbose
#undef mbottom
#undef alphabottom
#undef mtop
#undef alphatop
#undef mleft
#undef alphaleft
#undef mright
#undef alpharight
#undef W
#undef m
#undef alpha
#undef Qc
#undef R0
#undef option
#undef dimensionsAt
#undef majorAxisoffsetyh
#undef majorAxisoffsetxw
#undef minorAxisyh
#undef majorAxisyh
#undef minorAxisxw
#undef majorAxisxw
#undef loutyh
#undef linyh
#undef loutxw
#undef linxw
#undef yheight
#undef xwidth
#undef l
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component EndOfelement_1. */
  SIG_MESSAGE("EndOfelement_1 (Init)");

  /* Initializations for component SlitGuideEnd. */
  SIG_MESSAGE("SlitGuideEnd (Init)");
#define mccompcurname  SlitGuideEnd
#define mccompcurtype  Slit
#define mccompcurindex 13
#define xmin mccSlitGuideEnd_xmin
#define xmax mccSlitGuideEnd_xmax
#define ymin mccSlitGuideEnd_ymin
#define ymax mccSlitGuideEnd_ymax
#define radius mccSlitGuideEnd_radius
#define xwidth mccSlitGuideEnd_xwidth
#define yheight mccSlitGuideEnd_yheight
#line 50 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 14310 "./FullInstrument_v4.7_frontend.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component SlitMonitor. */
  SIG_MESSAGE("SlitMonitor (Init)");
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
#define nx mccSlitMonitor_nx
#define ny mccSlitMonitor_ny
#define filename mccSlitMonitor_filename
#define xmin mccSlitMonitor_xmin
#define xmax mccSlitMonitor_xmax
#define ymin mccSlitMonitor_ymin
#define ymax mccSlitMonitor_ymax
#define xwidth mccSlitMonitor_xwidth
#define yheight mccSlitMonitor_yheight
#define restore_neutron mccSlitMonitor_restore_neutron
#define nowritefile mccSlitMonitor_nowritefile
#line 69 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 14366 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component MonoFLat. */
  SIG_MESSAGE("MonoFLat (Init)");
#define mccompcurname  MonoFLat
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 15
#define mos_rms_y mccMonoFLat_mos_rms_y
#define mos_rms_z mccMonoFLat_mos_rms_z
#define mos_rms_max mccMonoFLat_mos_rms_max
#define mono_Q mccMonoFLat_mono_Q
#define SlabWidth mccMonoFLat_SlabWidth
#define SlabHeight mccMonoFLat_SlabHeight
#define rTable mccMonoFLat_rTable
#define tTable mccMonoFLat_tTable
#define row mccMonoFLat_row
#define col mccMonoFLat_col
#define tiltH mccMonoFLat_tiltH
#define tiltV mccMonoFLat_tiltV
#define reflect mccMonoFLat_reflect
#define transmit mccMonoFLat_transmit
#define zwidth mccMonoFLat_zwidth
#define yheight mccMonoFLat_yheight
#define gap mccMonoFLat_gap
#define NH mccMonoFLat_NH
#define NV mccMonoFLat_NV
#define mosaich mccMonoFLat_mosaich
#define mosaicv mccMonoFLat_mosaicv
#define r0 mccMonoFLat_r0
#define t0 mccMonoFLat_t0
#define Q mccMonoFLat_Q
#define RV mccMonoFLat_RV
#define RH mccMonoFLat_RH
#define DM mccMonoFLat_DM
#define mosaic mccMonoFLat_mosaic
#define width mccMonoFLat_width
#define height mccMonoFLat_height
#define verbose mccMonoFLat_verbose
#define order mccMonoFLat_order
#line 148 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  int i;

  if (mosaic != 0) {
    mos_rms_y = MIN2RAD*mosaic/sqrt(8*log(2));
    mos_rms_z = mos_rms_y; }
  else {
    mos_rms_y = MIN2RAD*mosaich/sqrt(8*log(2));
    mos_rms_z = MIN2RAD*mosaicv/sqrt(8*log(2)); }
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (mono_Q <= 0) { fprintf(stderr,"Monochromator_curved: %s: Error scattering vector Q = 0\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 <  0) { fprintf(stderr,"Monochromator_curved: %s: Error reflectivity r0 is negative\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 == 0) { fprintf(stderr,"Monochromator_curved: %s: Reflectivity r0 is null. Ignoring component.\n", NAME_CURRENT_COMP); }
  if (NH*NV == 0) { fprintf(stderr,"Monochromator_curved: %s: no slabs ??? (NH or NV=0)\n", NAME_CURRENT_COMP); exit(-1); }


  if (verbose && r0)
  {
    printf("Monochromator_curved: component %s Q=%.3g Angs-1 (DM=%.4g Angs)\n", NAME_CURRENT_COMP, mono_Q, 2*PI/mono_Q);
    if (NH*NV == 1) printf("            flat.\n");
    else
    { if (NH > 1)
      { printf("            horizontal: %i blades", (int)NH);
        if (RH != 0) printf(" focusing with RH=%.3g [m]", RH);
        printf("\n");
      }
      if (NV > 1)
      { printf("            vertical:   %i blades", (int)NV);
        if (RV != 0) printf(" focusing with RV=%.3g [m]", RV);
        printf("\n");
      }
    }
  }

  if (reflect != NULL && r0 && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Reflectivity data (k, R) from %s\n", NAME_CURRENT_COMP, reflect);
    Table_Read(&rTable, reflect, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&rTable);         /* rebin as evenly, increasing array */
    if (rTable.rows < 2) Table_Free(&rTable);
    if (verbose) Table_Info(rTable);
  } else rTable.data = NULL;

  if (transmit != NULL && strlen(transmit) && strcmp(transmit,"NULL") && strcmp(transmit,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Transmission data (k, T) from %s\n", NAME_CURRENT_COMP, transmit);
    Table_Read(&tTable, transmit, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&tTable);         /* rebin as evenly, increasing array */
    if (tTable.rows < 2) Table_Free(&tTable);
    if (verbose) Table_Info(tTable);
  } else tTable.data = NULL;

  if (width == 0) SlabWidth = zwidth;
  else SlabWidth = (width+gap)/NH - gap;
  if (height == 0) SlabHeight = yheight;
  else SlabHeight = (height+gap)/NV - gap;

  tiltH=malloc((int)(NH+1)*sizeof(double));
  tiltV=malloc((int)(NV+1)*sizeof(double));

  if (!tiltH) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NH=%g).\n", NAME_CURRENT_COMP, NH);
  else if (RH) { /* pre-compute tilts */
    for (i=0;i<=NH;i++)
    {
      tiltH[i]=asin((i-(NH+1)/2)*(SlabWidth+gap)/RH);
    }
  }
  if (!tiltV) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NV=%g).\n", NAME_CURRENT_COMP, NV);
  else if (RV) {
    for (i=0;i<=NV;i++)
    {
      tiltV[i]=-asin((i-(NV+1)/2)*(SlabHeight+gap)/RV);
    }
  }

}
#line 14503 "./FullInstrument_v4.7_frontend.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Mono. */
  SIG_MESSAGE("Mono (Init)");
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 16
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define SlabWidth mccMono_SlabWidth
#define SlabHeight mccMono_SlabHeight
#define rTable mccMono_rTable
#define tTable mccMono_tTable
#define row mccMono_row
#define col mccMono_col
#define tiltH mccMono_tiltH
#define tiltV mccMono_tiltV
#define reflect mccMono_reflect
#define transmit mccMono_transmit
#define zwidth mccMono_zwidth
#define yheight mccMono_yheight
#define gap mccMono_gap
#define NH mccMono_NH
#define NV mccMono_NV
#define mosaich mccMono_mosaich
#define mosaicv mccMono_mosaicv
#define r0 mccMono_r0
#define t0 mccMono_t0
#define Q mccMono_Q
#define RV mccMono_RV
#define RH mccMono_RH
#define DM mccMono_DM
#define mosaic mccMono_mosaic
#define width mccMono_width
#define height mccMono_height
#define verbose mccMono_verbose
#define order mccMono_order
#line 148 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  int i;

  if (mosaic != 0) {
    mos_rms_y = MIN2RAD*mosaic/sqrt(8*log(2));
    mos_rms_z = mos_rms_y; }
  else {
    mos_rms_y = MIN2RAD*mosaich/sqrt(8*log(2));
    mos_rms_z = MIN2RAD*mosaicv/sqrt(8*log(2)); }
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (mono_Q <= 0) { fprintf(stderr,"Monochromator_curved: %s: Error scattering vector Q = 0\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 <  0) { fprintf(stderr,"Monochromator_curved: %s: Error reflectivity r0 is negative\n", NAME_CURRENT_COMP); exit(-1); }
  if (r0 == 0) { fprintf(stderr,"Monochromator_curved: %s: Reflectivity r0 is null. Ignoring component.\n", NAME_CURRENT_COMP); }
  if (NH*NV == 0) { fprintf(stderr,"Monochromator_curved: %s: no slabs ??? (NH or NV=0)\n", NAME_CURRENT_COMP); exit(-1); }


  if (verbose && r0)
  {
    printf("Monochromator_curved: component %s Q=%.3g Angs-1 (DM=%.4g Angs)\n", NAME_CURRENT_COMP, mono_Q, 2*PI/mono_Q);
    if (NH*NV == 1) printf("            flat.\n");
    else
    { if (NH > 1)
      { printf("            horizontal: %i blades", (int)NH);
        if (RH != 0) printf(" focusing with RH=%.3g [m]", RH);
        printf("\n");
      }
      if (NV > 1)
      { printf("            vertical:   %i blades", (int)NV);
        if (RV != 0) printf(" focusing with RV=%.3g [m]", RV);
        printf("\n");
      }
    }
  }

  if (reflect != NULL && r0 && strlen(reflect) && strcmp(reflect,"NULL") && strcmp(reflect,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Reflectivity data (k, R) from %s\n", NAME_CURRENT_COMP, reflect);
    Table_Read(&rTable, reflect, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&rTable);         /* rebin as evenly, increasing array */
    if (rTable.rows < 2) Table_Free(&rTable);
    if (verbose) Table_Info(rTable);
  } else rTable.data = NULL;

  if (transmit != NULL && strlen(transmit) && strcmp(transmit,"NULL") && strcmp(transmit,"0"))
  {
    if (verbose) fprintf(stdout, "Monochromator_curved: %s: Transmission data (k, T) from %s\n", NAME_CURRENT_COMP, transmit);
    Table_Read(&tTable, transmit, 1); /* read 1st block data from file into rTable */
    Table_Rebin(&tTable);         /* rebin as evenly, increasing array */
    if (tTable.rows < 2) Table_Free(&tTable);
    if (verbose) Table_Info(tTable);
  } else tTable.data = NULL;

  if (width == 0) SlabWidth = zwidth;
  else SlabWidth = (width+gap)/NH - gap;
  if (height == 0) SlabHeight = yheight;
  else SlabHeight = (height+gap)/NV - gap;

  tiltH=malloc((int)(NH+1)*sizeof(double));
  tiltV=malloc((int)(NV+1)*sizeof(double));

  if (!tiltH) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NH=%g).\n", NAME_CURRENT_COMP, NH);
  else if (RH) { /* pre-compute tilts */
    for (i=0;i<=NH;i++)
    {
      tiltH[i]=asin((i-(NH+1)/2)*(SlabWidth+gap)/RH);
    }
  }
  if (!tiltV) printf("Monochromator_curved: %s: Warning: not enough memory to allocate tilts (NV=%g).\n", NAME_CURRENT_COMP, NV);
  else if (RV) {
    for (i=0;i<=NV;i++)
    {
      tiltV[i]=-asin((i-(NV+1)/2)*(SlabHeight+gap)/RV);
    }
  }

}
#line 14658 "./FullInstrument_v4.7_frontend.c"
#undef order
#undef verbose
#undef height
#undef width
#undef mosaic
#undef DM
#undef RH
#undef RV
#undef Q
#undef t0
#undef r0
#undef mosaicv
#undef mosaich
#undef NV
#undef NH
#undef gap
#undef yheight
#undef zwidth
#undef transmit
#undef reflect
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component A2. */
  SIG_MESSAGE("A2 (Init)");

  /* Initializations for component slit. */
  SIG_MESSAGE("slit (Init)");
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 18
#define xmin mccslit_xmin
#define xmax mccslit_xmax
#define ymin mccslit_ymin
#define ymax mccslit_ymax
#define radius mccslit_radius
#define xwidth mccslit_xwidth
#define yheight mccslit_yheight
#line 50 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 14730 "./FullInstrument_v4.7_frontend.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component SamplePos. */
  SIG_MESSAGE("SamplePos (Init)");

  /* Initializations for component psd_monitor_sample_pos. */
  SIG_MESSAGE("psd_monitor_sample_pos (Init)");
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
#define nx mccpsd_monitor_sample_pos_nx
#define ny mccpsd_monitor_sample_pos_ny
#define filename mccpsd_monitor_sample_pos_filename
#define xmin mccpsd_monitor_sample_pos_xmin
#define xmax mccpsd_monitor_sample_pos_xmax
#define ymin mccpsd_monitor_sample_pos_ymin
#define ymax mccpsd_monitor_sample_pos_ymax
#define xwidth mccpsd_monitor_sample_pos_xwidth
#define yheight mccpsd_monitor_sample_pos_yheight
#define restore_neutron mccpsd_monitor_sample_pos_restore_neutron
#define nowritefile mccpsd_monitor_sample_pos_nowritefile
#line 69 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 14789 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component e_monitor_sample_ps. */
  SIG_MESSAGE("e_monitor_sample_ps (Init)");
#define mccompcurname  e_monitor_sample_ps
#define mccompcurtype  E_monitor
#define mccompcurindex 21
#define nE mcce_monitor_sample_ps_nE
#define E_N mcce_monitor_sample_ps_E_N
#define E_p mcce_monitor_sample_ps_E_p
#define E_p2 mcce_monitor_sample_ps_E_p2
#define S_p mcce_monitor_sample_ps_S_p
#define S_pE mcce_monitor_sample_ps_S_pE
#define S_pE2 mcce_monitor_sample_ps_S_pE2
#define filename mcce_monitor_sample_ps_filename
#define xmin mcce_monitor_sample_ps_xmin
#define xmax mcce_monitor_sample_ps_xmax
#define ymin mcce_monitor_sample_ps_ymin
#define ymax mcce_monitor_sample_ps_ymax
#define xwidth mcce_monitor_sample_ps_xwidth
#define yheight mcce_monitor_sample_ps_yheight
#define Emin mcce_monitor_sample_ps_Emin
#define Emax mcce_monitor_sample_ps_Emax
#define restore_neutron mcce_monitor_sample_ps_restore_neutron
#define nowritefile mcce_monitor_sample_ps_nowritefile
#line 66 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
int i;

if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
    if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

    if ((xmin >= xmax) || (ymin >= ymax)) {
            printf("E_monitor: %s: Null detection area !\n"
                   "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
      exit(-1);
    }

    for (i=0; i<nE; i++)
    {
      E_N[i] = 0;
      E_p[i] = 0;
      E_p2[i] = 0;
    }
    S_p = S_pE = S_pE2 = 0;
}
#line 14853 "./FullInstrument_v4.7_frontend.c"
#undef nowritefile
#undef restore_neutron
#undef Emax
#undef Emin
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component Origin [1] */
  mccoordschange(mcposrOrigin, mcrotrOrigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Origin (without coords transformations) */
  mcJumpTrace_Origin:
  SIG_MESSAGE("Origin (Trace)");
  mcDEBUG_COMP("Origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompOrigin
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 70 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) mcsave(NULL);
  }
}
#line 15029 "./FullInstrument_v4.7_frontend.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompOrigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component source [2] */
  mccoordschange(mcposrsource, mcrotrsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source (without coords transformations) */
  mcJumpTrace_source:
  SIG_MESSAGE("source (Trace)");
  mcDEBUG_COMP("source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  source
#define mccompcurtype  Source_gen4
#define mccompcurindex 2
#define flux_file mccsource_flux_file
#define xdiv_file mccsource_xdiv_file
#define ydiv_file mccsource_ydiv_file
#define p_in mccsource_p_in
#define lambda0 mccsource_lambda0
#define lambda02 mccsource_lambda02
#define L2P mccsource_L2P
#define lambda0b mccsource_lambda0b
#define lambda02b mccsource_lambda02b
#define L2Pb mccsource_L2Pb
#define lambda0c mccsource_lambda0c
#define lambda02c mccsource_lambda02c
#define L2Pc mccsource_L2Pc
#define pTable mccsource_pTable
#define pTable_x mccsource_pTable_x
#define pTable_y mccsource_pTable_y
{   /* Declarations of source=Source_gen4() SETTING parameters. */
MCNUM radius = mccsource_radius;
MCNUM dist = mccsource_dist;
MCNUM xw = mccsource_xw;
MCNUM yh = mccsource_yh;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM Lambda0 = mccsource_Lambda0;
MCNUM dLambda = mccsource_dLambda;
MCNUM I1 = mccsource_I1;
MCNUM h = mccsource_h;
MCNUM w = mccsource_w;
MCNUM gaussian = mccsource_gaussian;
MCNUM verbose = mccsource_verbose;
MCNUM T1 = mccsource_T1;
MCNUM flux_file_perAA = mccsource_flux_file_perAA;
MCNUM flux_file_log = mccsource_flux_file_log;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM Emin = mccsource_Emin;
MCNUM Emax = mccsource_Emax;
MCNUM T2 = mccsource_T2;
MCNUM I2 = mccsource_I2;
MCNUM T3 = mccsource_T3;
MCNUM I3 = mccsource_I3;
MCNUM length = mccsource_length;
MCNUM phi_init = mccsource_phi_init;
MCNUM theta_init = mccsource_theta_init;
MCNUM HEtailA = mccsource_HEtailA;
MCNUM HEtailL0 = mccsource_HEtailL0;
#line 451 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
{
  double theta0,phi0,theta1,phi1,chi,theta,phi,v,r, lambda;
  double tan_h, tan_v, Maxwell, lambda2, lambda5;

  if (verbose >= 0)
  {

    z=0;

    if ((h == 0) || (w == 0))
    {
      chi=2*PI*rand01();                          /* Choose point on source */
      r=sqrt(rand01())*radius;                    /* with uniform distribution. */
      x=r*cos(chi);
      y=r*sin(chi);
    }
    else
    {
      x = w*randpm1()/2;   /* select point on source (uniform) */
      y = h*randpm1()/2;
    }
    if (length != 0)
      z = length*randpm1()/2;

    if (dist == 0)
    {
      theta0 =  DEG2RAD*xw/2;
      phi0   =  DEG2RAD*yh/2;
      theta1 = -DEG2RAD*xw/2;
      phi1   = -DEG2RAD*yh/2;
    }
    else
    {
      theta0= -atan((x-xw/2.0)/dist)+theta_init/180*3.1415;              /* Angles to aim at target */
      phi0  = -atan((y-yh/2.0)/dist)+phi_init/180*3.1415;
      theta1= -atan((x+xw/2.0)/dist)+theta_init/180*3.1415;
      phi1  = -atan((y+yh/2.0)/dist)+phi_init/180*3.1415;
    }

    /* shot towards target : flat distribution */

    if (gaussian)
    {
      theta= theta0+(theta1- theta0)*(randnorm()*FWHM2RMS+0.5);
      phi  = phi0  +(phi1  - phi0)  *(randnorm()*FWHM2RMS+0.5);
    }
    else
    {
      theta= theta0+(theta1- theta0)*rand01();
      phi  = phi0  +(phi1  - phi0)  *rand01();
    }
    /* Assume linear distribution */
    lambda = Lambda0+dLambda*randpm1();
    if (lambda <= 0) ABSORB;

    v = K2V*(2*PI/lambda);
    p = p_in;
    if (!flux_file  || !strlen(flux_file))
      p *=  2 * fabs((theta1 - theta0)*sin((phi1 - phi0)/2)); /* solid angle */
    p *= cos(phi)*cos(theta);
    if (flux_file && strlen(flux_file) > 0)
    {
       double w=Table_Value(pTable, lambda, 1);
       if (flux_file_log) w=exp(w);
       p *= w;
    }
    else if (T1 > 0 && I1 > 0)
    {
      lambda2 = lambda*lambda;
      lambda5 = lambda2*lambda2*lambda;
      Maxwell = I1 * L2P/lambda5 * exp(-lambda02/lambda2);  /* 1/AA */

      if ((T2 > 0) && (I2 > 0))
      {
        Maxwell += I2 * L2Pb/lambda5 * exp(-lambda02b/lambda2);
      }
      if ((T3 > 0) && (I3 > 0))
      {
        Maxwell += I3 * L2Pc/lambda5 * exp(-lambda02c/lambda2);
      }
      if (HEtailA>0)
      {
        Maxwell+=HEtailA/(lambda-HEtailL0)/(lambda-HEtailL0);
      }
      p *= Maxwell;
    }

    /* Perform the correct treatment - no small angle approx. here! */
    tan_h = tan(theta);
    tan_v = tan(phi);
    vz = v / sqrt(1 + tan_v*tan_v + tan_h*tan_h);
    vy = tan_v * vz;
    vx = tan_h * vz;
    /* optional x-xdiv and y-ydiv weightening: position=along columns, div=along rows */
    if (xdiv_file && strlen(xdiv_file) > 0 && pTable_xsum > 0) {
      double i,j;
      j = (x-            pTable_xmin) /(pTable_xmax -pTable_xmin) *pTable_x.columns;
      i = (theta*RAD2DEG-pTable_dxmin)/(pTable_dxmax-pTable_dxmin)*pTable_x.rows;
      r = Table_Value2d(pTable_x, i,j); /* row, column */
      p *= r/pTable_xsum;
    }
    if (ydiv_file && strlen(ydiv_file) > 0 && pTable_ysum > 0) {
      double i,j;
      j = (y-            pTable_ymin) /(pTable_ymax -pTable_ymin) *pTable_y.columns;
      i = (phi*RAD2DEG-  pTable_dymin)/(pTable_dymax-pTable_dymin)*pTable_y.rows;
      r = Table_Value2d(pTable_y, i,j);
      p *= r/pTable_ysum;
    }
    SCATTER;
  }
}
#line 15296 "./FullInstrument_v4.7_frontend.c"
}   /* End of source=Source_gen4() SETTING parameter declarations. */
#undef pTable_y
#undef pTable_x
#undef pTable
#undef L2Pc
#undef lambda02c
#undef lambda0c
#undef L2Pb
#undef lambda02b
#undef lambda0b
#undef L2P
#undef lambda02
#undef lambda0
#undef p_in
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component StartOfGuide [3] */
  mccoordschange(mcposrStartOfGuide, mcrotrStartOfGuide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component StartOfGuide (without coords transformations) */
  mcJumpTrace_StartOfGuide:
  SIG_MESSAGE("StartOfGuide (Trace)");
  mcDEBUG_COMP("StartOfGuide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompStartOfGuide
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  StartOfGuide
#define mccompcurtype  Arm
#define mccompcurindex 3
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompStartOfGuide:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component elliptical_guide_gravity5 [4] */
  mccoordschange(mcposrelliptical_guide_gravity5, mcrotrelliptical_guide_gravity5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component elliptical_guide_gravity5 (without coords transformations) */
  mcJumpTrace_elliptical_guide_gravity5:
  SIG_MESSAGE("elliptical_guide_gravity5 (Trace)");
  mcDEBUG_COMP("elliptical_guide_gravity5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompelliptical_guide_gravity5
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  elliptical_guide_gravity5
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 4
#define mvaluesright mccelliptical_guide_gravity5_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity5_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity5_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity5_mvaluesbottom
#define seglength mccelliptical_guide_gravity5_seglength
#define guideInfo mccelliptical_guide_gravity5_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity5_latestParticleCollision
#define Gx mccelliptical_guide_gravity5_Gx
#define Gy mccelliptical_guide_gravity5_Gy
#define Gz mccelliptical_guide_gravity5_Gz
#define Gx0 mccelliptical_guide_gravity5_Gx0
#define Gy0 mccelliptical_guide_gravity5_Gy0
#define Gz0 mccelliptical_guide_gravity5_Gz0
#define Circ mccelliptical_guide_gravity5_Circ
#define dynamicalSegLength mccelliptical_guide_gravity5_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity5=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity5_l;
MCNUM xwidth = mccelliptical_guide_gravity5_xwidth;
MCNUM yheight = mccelliptical_guide_gravity5_yheight;
MCNUM linxw = mccelliptical_guide_gravity5_linxw;
MCNUM loutxw = mccelliptical_guide_gravity5_loutxw;
MCNUM linyh = mccelliptical_guide_gravity5_linyh;
MCNUM loutyh = mccelliptical_guide_gravity5_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity5_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity5_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity5_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity5_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity5_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity5_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity5_dimensionsAt;
char* option = mccelliptical_guide_gravity5_option;
MCNUM R0 = mccelliptical_guide_gravity5_R0;
MCNUM Qc = mccelliptical_guide_gravity5_Qc;
MCNUM alpha = mccelliptical_guide_gravity5_alpha;
MCNUM m = mccelliptical_guide_gravity5_m;
MCNUM W = mccelliptical_guide_gravity5_W;
MCNUM alpharight = mccelliptical_guide_gravity5_alpharight;
MCNUM mright = mccelliptical_guide_gravity5_mright;
MCNUM alphaleft = mccelliptical_guide_gravity5_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity5_mleft;
MCNUM alphatop = mccelliptical_guide_gravity5_alphatop;
MCNUM mtop = mccelliptical_guide_gravity5_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity5_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity5_mbottom;
char* verbose = mccelliptical_guide_gravity5_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity5_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity5_curvature;
#line 1420 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
  PROP_Z0;
  SCATTER;
  
  double Gloc;
  if (curvature) {
    Gloc=(vx*vx+vy*vy+vz*vz)/curvature;
  } else {
    Gloc=0;
  }
  

	if( !guideInfo.EnclosingBoxOn )
		if( fabs(x) > guideInfo.entranceHorizontalWidth/2.0 
		    || fabs(y) > guideInfo.entranceVerticalWidth/2.0 )
		  ABSORB;


	int bounces = 0;
	for(bounces = 0; bounces <= 1000; bounces++){
	
	  Gx=Gx0; Gy=Gy0; Gz=Gz0;
	  if (curvature) {
	    // Add velocity-dependent, location-dependent approximation to centripetal force for curvature...
	    Gx=Gx0+Gloc*cos(2*PI*z/Circ);Gz=Gz0+Gloc*sin(2*PI*z/Circ);
	  }

		// Find the next intersection between the guide and the neutron.
		int boolean = guide_elliptical_handleGuideIntersection(
							x,y,z,vx,vy,vz,Gx,Gy,Gz,
							&guideInfo,&latestParticleCollision);

		double timeToCollision = 
				latestParticleCollision.delta_time_to_next_collision;

		// Handle special cases.
		if(boolean == 0) ABSORB;
		if(timeToCollision < 1e-15) ABSORB;

		// If the neutron reach the end of the guide, when move 
		// the neutron to the end of guide and leave this component.
		if( z+vz*timeToCollision+0.5*Gz*timeToCollision*timeToCollision 
			>= guideInfo.Length )
		{
			double timeToExit = 0;
			solve_2nd_order(
					&timeToExit,NULL,
					-0.5*Gz,-vz,guideInfo.Length-z-1e-9);
			PROP_GRAV_DT(timeToExit,Gx,Gy,Gz);
			SCATTER;
			break;
		}

		// Move the neutron and handle the reflection.
		PROP_GRAV_DT(timeToCollision,Gx,Gy,Gz);
		if(latestParticleCollision.collisionType == Absorb){ ABSORB; }
		if(latestParticleCollision.collisionType == Reflex){
			p *= guide_elliptical_handleReflection(x,y,z,&vx,&vy,&vz,
						&guideInfo,&latestParticleCollision); 
			SCATTER;
			if(p == 0) ABSORB; 
		}
	}
	
	if( fabs(x) > guideInfo.exitHorizontalWidth/2 
	 || fabs(y) > guideInfo.exitVerticalWidth/2 )
		ABSORB;
	
}
#line 15637 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity5=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompelliptical_guide_gravity5:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component EndOfelement_5 [5] */
  mccoordschange(mcposrEndOfelement_5, mcrotrEndOfelement_5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EndOfelement_5 (without coords transformations) */
  mcJumpTrace_EndOfelement_5:
  SIG_MESSAGE("EndOfelement_5 (Trace)");
  mcDEBUG_COMP("EndOfelement_5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompEndOfelement_5
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  EndOfelement_5
#define mccompcurtype  Arm
#define mccompcurindex 5
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompEndOfelement_5:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component EndOfelement_4 [6] */
  mccoordschange(mcposrEndOfelement_4, mcrotrEndOfelement_4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EndOfelement_4 (without coords transformations) */
  mcJumpTrace_EndOfelement_4:
  SIG_MESSAGE("EndOfelement_4 (Trace)");
  mcDEBUG_COMP("EndOfelement_4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompEndOfelement_4
  STORE_NEUTRON(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  EndOfelement_4
#define mccompcurtype  Arm
#define mccompcurindex 6
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompEndOfelement_4:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component elliptical_guide_gravity3 [7] */
  mccoordschange(mcposrelliptical_guide_gravity3, mcrotrelliptical_guide_gravity3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component elliptical_guide_gravity3 (without coords transformations) */
  mcJumpTrace_elliptical_guide_gravity3:
  SIG_MESSAGE("elliptical_guide_gravity3 (Trace)");
  mcDEBUG_COMP("elliptical_guide_gravity3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompelliptical_guide_gravity3
  STORE_NEUTRON(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  elliptical_guide_gravity3
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 7
#define mvaluesright mccelliptical_guide_gravity3_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity3_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity3_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity3_mvaluesbottom
#define seglength mccelliptical_guide_gravity3_seglength
#define guideInfo mccelliptical_guide_gravity3_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity3_latestParticleCollision
#define Gx mccelliptical_guide_gravity3_Gx
#define Gy mccelliptical_guide_gravity3_Gy
#define Gz mccelliptical_guide_gravity3_Gz
#define Gx0 mccelliptical_guide_gravity3_Gx0
#define Gy0 mccelliptical_guide_gravity3_Gy0
#define Gz0 mccelliptical_guide_gravity3_Gz0
#define Circ mccelliptical_guide_gravity3_Circ
#define dynamicalSegLength mccelliptical_guide_gravity3_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity3=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity3_l;
MCNUM xwidth = mccelliptical_guide_gravity3_xwidth;
MCNUM yheight = mccelliptical_guide_gravity3_yheight;
MCNUM linxw = mccelliptical_guide_gravity3_linxw;
MCNUM loutxw = mccelliptical_guide_gravity3_loutxw;
MCNUM linyh = mccelliptical_guide_gravity3_linyh;
MCNUM loutyh = mccelliptical_guide_gravity3_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity3_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity3_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity3_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity3_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity3_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity3_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity3_dimensionsAt;
char* option = mccelliptical_guide_gravity3_option;
MCNUM R0 = mccelliptical_guide_gravity3_R0;
MCNUM Qc = mccelliptical_guide_gravity3_Qc;
MCNUM alpha = mccelliptical_guide_gravity3_alpha;
MCNUM m = mccelliptical_guide_gravity3_m;
MCNUM W = mccelliptical_guide_gravity3_W;
MCNUM alpharight = mccelliptical_guide_gravity3_alpharight;
MCNUM mright = mccelliptical_guide_gravity3_mright;
MCNUM alphaleft = mccelliptical_guide_gravity3_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity3_mleft;
MCNUM alphatop = mccelliptical_guide_gravity3_alphatop;
MCNUM mtop = mccelliptical_guide_gravity3_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity3_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity3_mbottom;
char* verbose = mccelliptical_guide_gravity3_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity3_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity3_curvature;
#line 1420 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
  PROP_Z0;
  SCATTER;
  
  double Gloc;
  if (curvature) {
    Gloc=(vx*vx+vy*vy+vz*vz)/curvature;
  } else {
    Gloc=0;
  }
  

	if( !guideInfo.EnclosingBoxOn )
		if( fabs(x) > guideInfo.entranceHorizontalWidth/2.0 
		    || fabs(y) > guideInfo.entranceVerticalWidth/2.0 )
		  ABSORB;


	int bounces = 0;
	for(bounces = 0; bounces <= 1000; bounces++){
	
	  Gx=Gx0; Gy=Gy0; Gz=Gz0;
	  if (curvature) {
	    // Add velocity-dependent, location-dependent approximation to centripetal force for curvature...
	    Gx=Gx0+Gloc*cos(2*PI*z/Circ);Gz=Gz0+Gloc*sin(2*PI*z/Circ);
	  }

		// Find the next intersection between the guide and the neutron.
		int boolean = guide_elliptical_handleGuideIntersection(
							x,y,z,vx,vy,vz,Gx,Gy,Gz,
							&guideInfo,&latestParticleCollision);

		double timeToCollision = 
				latestParticleCollision.delta_time_to_next_collision;

		// Handle special cases.
		if(boolean == 0) ABSORB;
		if(timeToCollision < 1e-15) ABSORB;

		// If the neutron reach the end of the guide, when move 
		// the neutron to the end of guide and leave this component.
		if( z+vz*timeToCollision+0.5*Gz*timeToCollision*timeToCollision 
			>= guideInfo.Length )
		{
			double timeToExit = 0;
			solve_2nd_order(
					&timeToExit,NULL,
					-0.5*Gz,-vz,guideInfo.Length-z-1e-9);
			PROP_GRAV_DT(timeToExit,Gx,Gy,Gz);
			SCATTER;
			break;
		}

		// Move the neutron and handle the reflection.
		PROP_GRAV_DT(timeToCollision,Gx,Gy,Gz);
		if(latestParticleCollision.collisionType == Absorb){ ABSORB; }
		if(latestParticleCollision.collisionType == Reflex){
			p *= guide_elliptical_handleReflection(x,y,z,&vx,&vy,&vz,
						&guideInfo,&latestParticleCollision); 
			SCATTER;
			if(p == 0) ABSORB; 
		}
	}
	
	if( fabs(x) > guideInfo.exitHorizontalWidth/2 
	 || fabs(y) > guideInfo.exitVerticalWidth/2 )
		ABSORB;
	
}
#line 16080 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity3=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompelliptical_guide_gravity3:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component EndOfelement_3 [8] */
  mccoordschange(mcposrEndOfelement_3, mcrotrEndOfelement_3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EndOfelement_3 (without coords transformations) */
  mcJumpTrace_EndOfelement_3:
  SIG_MESSAGE("EndOfelement_3 (Trace)");
  mcDEBUG_COMP("EndOfelement_3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompEndOfelement_3
  STORE_NEUTRON(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  EndOfelement_3
#define mccompcurtype  Arm
#define mccompcurindex 8
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompEndOfelement_3:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Curved_guide_2 [9] */
  mccoordschange(mcposrCurved_guide_2, mcrotrCurved_guide_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Curved_guide_2 (without coords transformations) */
  mcJumpTrace_Curved_guide_2:
  SIG_MESSAGE("Curved_guide_2 (Trace)");
  mcDEBUG_COMP("Curved_guide_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompCurved_guide_2
  STORE_NEUTRON(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  Curved_guide_2
#define mccompcurtype  Bender
#define mccompcurindex 9
#define bk mccCurved_guide_2_bk
#define mWin mccCurved_guide_2_mWin
{   /* Declarations of Curved_guide_2=Bender() SETTING parameters. */
MCNUM w = mccCurved_guide_2_w;
MCNUM h = mccCurved_guide_2_h;
MCNUM r = mccCurved_guide_2_r;
MCNUM Win = mccCurved_guide_2_Win;
MCNUM k = mccCurved_guide_2_k;
MCNUM d = mccCurved_guide_2_d;
MCNUM l = mccCurved_guide_2_l;
MCNUM R0a = mccCurved_guide_2_R0a;
MCNUM Qca = mccCurved_guide_2_Qca;
MCNUM alphaa = mccCurved_guide_2_alphaa;
MCNUM ma = mccCurved_guide_2_ma;
MCNUM Wa = mccCurved_guide_2_Wa;
MCNUM R0i = mccCurved_guide_2_R0i;
MCNUM Qci = mccCurved_guide_2_Qci;
MCNUM alphai = mccCurved_guide_2_alphai;
MCNUM mi = mccCurved_guide_2_mi;
MCNUM Wi = mccCurved_guide_2_Wi;
MCNUM R0s = mccCurved_guide_2_R0s;
MCNUM Qcs = mccCurved_guide_2_Qcs;
MCNUM alphas = mccCurved_guide_2_alphas;
MCNUM ms = mccCurved_guide_2_ms;
MCNUM Ws = mccCurved_guide_2_Ws;
#line 133 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
    int i,num,numa,numi;
    double dru,ab,dab,R,Q,Ta,vpl;
    double einmWin,ausmWin,zykmWin,aeumWin,innmWin,ref,innref,aeuref;
    double einzei,auszei,zykzei;

    /* does the neutron hit the bender at the entrance? */
    PROP_Z0;
    if ((fabs(x)<w/2) && (fabs(y)<h/2))
    {
      /*** reflections in the XZ-plane ***/

      /* distance between neutron and concave side of the channel at the entrance */
      dru=floor((w/2-x)/bk)*bk;
      ab=w/2.0-x-dru;

      /* radius of the channel */
      R=r-dru;

      /* does the neutron hit the partition at the entrance? */
      if (ab<bk-d)
      {
        double aeu[] = {R0a, Qca, alphaa, ma, Wa};
        /* velocity in the XZ-plane */
        vpl=sqrt(vx*vx+vz*vz);

        /* divergence of the neutron at the entrance */
        einmWin=atan(vx/vz);

        /* maximal distance between neutron and concave side of the channel */
        dab=R-cos(einmWin)*(R-ab);

        /* reflection angle at the concave side */
        aeumWin=acos((R-dab)/R);

        /* reflection coefficient at the concave side */
        Q=2.0*V2K*vpl*sin(aeumWin);
        StdReflecFunc(Q, aeu, &aeuref);

        /* does the neutron hit the convex side of the channel? */
        innmWin=0.0;
        innref=1.0;
        if (dab>bk-d)
        {
           double inn[] = {R0i, Qci, alphai, mi, Wi};
           /* reflection coefficient at the convex side */
           innmWin=acos((R-dab)/(R-bk+d));
           Q=2.0*V2K*vpl*sin(innmWin);
           StdReflecFunc(Q, inn, &innref);
        }

        /* divergence of the neutron at the exit */
        zykmWin=2.0*(aeumWin-innmWin);
        ausmWin=fmod(mWin+einmWin+aeumWin-innmWin
          *(1.0+SIGN(einmWin)),zykmWin)-zykmWin/2.0;
        ausmWin+=innmWin*SIGN(ausmWin);

        /* number of reflections at the concave side */
        numa=(mWin+einmWin+aeumWin-innmWin*(1.0+SIGN(einmWin)))/zykmWin;

        /* number of reflections at the convex side */
        numi=numa;
        if (ausmWin*einmWin<0)
        {
           if (ausmWin-einmWin>0)
              numi++;
           else
              numi--;
        }

        /* is the reflection coefficient too small? */
        if (((numa>0) && (aeuref<=0)) || ((numi>0) && (innref<=0)))
           ABSORB;

        /* calculation of the neutron probability weight p */
        for (i=1;i<=numa;i++)
            p*=aeuref;
        for (i=1;i<=numi;i++)
            p*=innref;

        /* time to cross the bender */
        Ta=(2*numa*(tan(aeumWin)-tan(innmWin))
          +tan(ausmWin)-tan(einmWin)
          -tan(innmWin)*(SIGN(ausmWin)-SIGN(einmWin)))
          *(R-dab)/vpl;
        t+=Ta;

        /* distance between neutron and concave side of channel at the exit */
        ab=R-(R-dab)/cos(ausmWin);

        /* calculation of the exit coordinates in the XZ-plane */
        x=w/2.0-ab-dru;
        z=r*mWin;
        vx=sin(ausmWin)*vpl;
        vz=cos(ausmWin)*vpl;

        /*** reflections at top and bottom side (Y axis) ***/

        if (vy!=0.0)
        {
          double s[] = {R0s, Qcs, alphas, ms, Ws};
          /* reflection coefficent at the top and bottom side */
          Q=2.0*V2K*fabs(vy);
          StdReflecFunc(Q, s, &ref);

          /* number of reflections at top and bottom */
          einzei=h/2.0/fabs(vy)+y/vy;
          zykzei=h/fabs(vy);
          num=(Ta+einzei)/zykzei;

          /* time between the last reflection and the exit */
          auszei=fmod(Ta+einzei,zykzei);

          /* is the reflection coefficient too small? */
          if ((num>0) && (ref<=0))
             ABSORB;

          /* calculation of the probability weight p */
          for (i=1;i<=num;i++) {
               p*=ref;
               vy*=-1.0; }

          /* calculation of the exit coordinate */
          y=auszei*vy-vy*h/fabs(vy)/2.0;
        } /* if (vy!=0.0) */
        SCATTER;
      } /* if (dab>bk-d)  */
      else
        ABSORB; /* hit separating walls */
    }
    else /* if ((fabs(x)<w/2) && (fabs(y)<h/2))   */
      ABSORB; /* miss entry window */

}
#line 16463 "./FullInstrument_v4.7_frontend.c"
}   /* End of Curved_guide_2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompCurved_guide_2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component EndOfelement_2 [10] */
  mccoordschange(mcposrEndOfelement_2, mcrotrEndOfelement_2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EndOfelement_2 (without coords transformations) */
  mcJumpTrace_EndOfelement_2:
  SIG_MESSAGE("EndOfelement_2 (Trace)");
  mcDEBUG_COMP("EndOfelement_2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompEndOfelement_2
  STORE_NEUTRON(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  EndOfelement_2
#define mccompcurtype  Arm
#define mccompcurindex 10
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompEndOfelement_2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component elliptical_guide_gravity1 [11] */
  mccoordschange(mcposrelliptical_guide_gravity1, mcrotrelliptical_guide_gravity1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component elliptical_guide_gravity1 (without coords transformations) */
  mcJumpTrace_elliptical_guide_gravity1:
  SIG_MESSAGE("elliptical_guide_gravity1 (Trace)");
  mcDEBUG_COMP("elliptical_guide_gravity1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompelliptical_guide_gravity1
  STORE_NEUTRON(11,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[11]++;
  mcPCounter[11] += p;
  mcP2Counter[11] += p*p;
#define mccompcurname  elliptical_guide_gravity1
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 11
#define mvaluesright mccelliptical_guide_gravity1_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity1_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity1_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity1_mvaluesbottom
#define seglength mccelliptical_guide_gravity1_seglength
#define guideInfo mccelliptical_guide_gravity1_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity1_latestParticleCollision
#define Gx mccelliptical_guide_gravity1_Gx
#define Gy mccelliptical_guide_gravity1_Gy
#define Gz mccelliptical_guide_gravity1_Gz
#define Gx0 mccelliptical_guide_gravity1_Gx0
#define Gy0 mccelliptical_guide_gravity1_Gy0
#define Gz0 mccelliptical_guide_gravity1_Gz0
#define Circ mccelliptical_guide_gravity1_Circ
#define dynamicalSegLength mccelliptical_guide_gravity1_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity1=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity1_l;
MCNUM xwidth = mccelliptical_guide_gravity1_xwidth;
MCNUM yheight = mccelliptical_guide_gravity1_yheight;
MCNUM linxw = mccelliptical_guide_gravity1_linxw;
MCNUM loutxw = mccelliptical_guide_gravity1_loutxw;
MCNUM linyh = mccelliptical_guide_gravity1_linyh;
MCNUM loutyh = mccelliptical_guide_gravity1_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity1_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity1_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity1_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity1_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity1_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity1_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity1_dimensionsAt;
char* option = mccelliptical_guide_gravity1_option;
MCNUM R0 = mccelliptical_guide_gravity1_R0;
MCNUM Qc = mccelliptical_guide_gravity1_Qc;
MCNUM alpha = mccelliptical_guide_gravity1_alpha;
MCNUM m = mccelliptical_guide_gravity1_m;
MCNUM W = mccelliptical_guide_gravity1_W;
MCNUM alpharight = mccelliptical_guide_gravity1_alpharight;
MCNUM mright = mccelliptical_guide_gravity1_mright;
MCNUM alphaleft = mccelliptical_guide_gravity1_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity1_mleft;
MCNUM alphatop = mccelliptical_guide_gravity1_alphatop;
MCNUM mtop = mccelliptical_guide_gravity1_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity1_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity1_mbottom;
char* verbose = mccelliptical_guide_gravity1_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity1_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity1_curvature;
#line 1420 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{
  PROP_Z0;
  SCATTER;
  
  double Gloc;
  if (curvature) {
    Gloc=(vx*vx+vy*vy+vz*vz)/curvature;
  } else {
    Gloc=0;
  }
  

	if( !guideInfo.EnclosingBoxOn )
		if( fabs(x) > guideInfo.entranceHorizontalWidth/2.0 
		    || fabs(y) > guideInfo.entranceVerticalWidth/2.0 )
		  ABSORB;


	int bounces = 0;
	for(bounces = 0; bounces <= 1000; bounces++){
	
	  Gx=Gx0; Gy=Gy0; Gz=Gz0;
	  if (curvature) {
	    // Add velocity-dependent, location-dependent approximation to centripetal force for curvature...
	    Gx=Gx0+Gloc*cos(2*PI*z/Circ);Gz=Gz0+Gloc*sin(2*PI*z/Circ);
	  }

		// Find the next intersection between the guide and the neutron.
		int boolean = guide_elliptical_handleGuideIntersection(
							x,y,z,vx,vy,vz,Gx,Gy,Gz,
							&guideInfo,&latestParticleCollision);

		double timeToCollision = 
				latestParticleCollision.delta_time_to_next_collision;

		// Handle special cases.
		if(boolean == 0) ABSORB;
		if(timeToCollision < 1e-15) ABSORB;

		// If the neutron reach the end of the guide, when move 
		// the neutron to the end of guide and leave this component.
		if( z+vz*timeToCollision+0.5*Gz*timeToCollision*timeToCollision 
			>= guideInfo.Length )
		{
			double timeToExit = 0;
			solve_2nd_order(
					&timeToExit,NULL,
					-0.5*Gz,-vz,guideInfo.Length-z-1e-9);
			PROP_GRAV_DT(timeToExit,Gx,Gy,Gz);
			SCATTER;
			break;
		}

		// Move the neutron and handle the reflection.
		PROP_GRAV_DT(timeToCollision,Gx,Gy,Gz);
		if(latestParticleCollision.collisionType == Absorb){ ABSORB; }
		if(latestParticleCollision.collisionType == Reflex){
			p *= guide_elliptical_handleReflection(x,y,z,&vx,&vy,&vz,
						&guideInfo,&latestParticleCollision); 
			SCATTER;
			if(p == 0) ABSORB; 
		}
	}
	
	if( fabs(x) > guideInfo.exitHorizontalWidth/2 
	 || fabs(y) > guideInfo.exitVerticalWidth/2 )
		ABSORB;
	
}
#line 16790 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity1=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompelliptical_guide_gravity1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(11,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component EndOfelement_1 [12] */
  mccoordschange(mcposrEndOfelement_1, mcrotrEndOfelement_1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component EndOfelement_1 (without coords transformations) */
  mcJumpTrace_EndOfelement_1:
  SIG_MESSAGE("EndOfelement_1 (Trace)");
  mcDEBUG_COMP("EndOfelement_1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompEndOfelement_1
  STORE_NEUTRON(12,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[12]++;
  mcPCounter[12] += p;
  mcP2Counter[12] += p*p;
#define mccompcurname  EndOfelement_1
#define mccompcurtype  Arm
#define mccompcurindex 12
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompEndOfelement_1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(12,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component SlitGuideEnd [13] */
  mccoordschange(mcposrSlitGuideEnd, mcrotrSlitGuideEnd,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SlitGuideEnd (without coords transformations) */
  mcJumpTrace_SlitGuideEnd:
  SIG_MESSAGE("SlitGuideEnd (Trace)");
  mcDEBUG_COMP("SlitGuideEnd")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSlitGuideEnd
  STORE_NEUTRON(13,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[13]++;
  mcPCounter[13] += p;
  mcP2Counter[13] += p*p;
#define mccompcurname  SlitGuideEnd
#define mccompcurtype  Slit
#define mccompcurindex 13
{   /* Declarations of SlitGuideEnd=Slit() SETTING parameters. */
MCNUM xmin = mccSlitGuideEnd_xmin;
MCNUM xmax = mccSlitGuideEnd_xmax;
MCNUM ymin = mccSlitGuideEnd_ymin;
MCNUM ymax = mccSlitGuideEnd_ymax;
MCNUM radius = mccSlitGuideEnd_radius;
MCNUM xwidth = mccSlitGuideEnd_xwidth;
MCNUM yheight = mccSlitGuideEnd_yheight;
#line 71 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 17032 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitGuideEnd=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSlitGuideEnd:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(13,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component SlitMonitor [14] */
  mccoordschange(mcposrSlitMonitor, mcrotrSlitMonitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SlitMonitor (without coords transformations) */
  mcJumpTrace_SlitMonitor:
  SIG_MESSAGE("SlitMonitor (Trace)");
  mcDEBUG_COMP("SlitMonitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSlitMonitor
  STORE_NEUTRON(14,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[14]++;
  mcPCounter[14] += p;
  mcP2Counter[14] += p*p;
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
{   /* Declarations of SlitMonitor=PSD_monitor() SETTING parameters. */
int nx = mccSlitMonitor_nx;
int ny = mccSlitMonitor_ny;
char* filename = mccSlitMonitor_filename;
MCNUM xmin = mccSlitMonitor_xmin;
MCNUM xmax = mccSlitMonitor_xmax;
MCNUM ymin = mccSlitMonitor_ymin;
MCNUM ymax = mccSlitMonitor_ymax;
MCNUM xwidth = mccSlitMonitor_xwidth;
MCNUM yheight = mccSlitMonitor_yheight;
MCNUM restore_neutron = mccSlitMonitor_restore_neutron;
int nowritefile = mccSlitMonitor_nowritefile;
#line 95 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 17167 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitMonitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSlitMonitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(14,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component MonoFLat [15] */
  mccoordschange(mcposrMonoFLat, mcrotrMonoFLat,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component MonoFLat (without coords transformations) */
  mcJumpTrace_MonoFLat:
  SIG_MESSAGE("MonoFLat (Trace)");
  mcDEBUG_COMP("MonoFLat")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMonoFLat
  STORE_NEUTRON(15,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[15]++;
  mcPCounter[15] += p;
  mcP2Counter[15] += p*p;
#define mccompcurname  MonoFLat
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 15
#define mos_rms_y mccMonoFLat_mos_rms_y
#define mos_rms_z mccMonoFLat_mos_rms_z
#define mos_rms_max mccMonoFLat_mos_rms_max
#define mono_Q mccMonoFLat_mono_Q
#define SlabWidth mccMonoFLat_SlabWidth
#define SlabHeight mccMonoFLat_SlabHeight
#define rTable mccMonoFLat_rTable
#define tTable mccMonoFLat_tTable
#define row mccMonoFLat_row
#define col mccMonoFLat_col
#define tiltH mccMonoFLat_tiltH
#define tiltV mccMonoFLat_tiltV
{   /* Declarations of MonoFLat=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonoFLat_reflect;
char* transmit = mccMonoFLat_transmit;
MCNUM zwidth = mccMonoFLat_zwidth;
MCNUM yheight = mccMonoFLat_yheight;
MCNUM gap = mccMonoFLat_gap;
MCNUM NH = mccMonoFLat_NH;
MCNUM NV = mccMonoFLat_NV;
MCNUM mosaich = mccMonoFLat_mosaich;
MCNUM mosaicv = mccMonoFLat_mosaicv;
MCNUM r0 = mccMonoFLat_r0;
MCNUM t0 = mccMonoFLat_t0;
MCNUM Q = mccMonoFLat_Q;
MCNUM RV = mccMonoFLat_RV;
MCNUM RH = mccMonoFLat_RH;
MCNUM DM = mccMonoFLat_DM;
MCNUM mosaic = mccMonoFLat_mosaic;
MCNUM width = mccMonoFLat_width;
MCNUM height = mccMonoFLat_height;
MCNUM verbose = mccMonoFLat_verbose;
MCNUM order = mccMonoFLat_order;
/* 'MonoFLat=Monochromator_curved()' component instance has conditional execution */
if (( mcipold_mono == 1 ))

#line 230 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  double dt;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0 && r0)
  {                             /* Moving towards crystal? */
    double zmin,zmax, ymin,ymax;

    /* Propagate to crystal plane */
    PROP_DT(dt);    /* now in the vertical plane of monochromator */

    zmax = ((NH*(SlabWidth+gap))-gap)/2;
    zmin = -zmax;
    ymax = ((NV*(SlabHeight+gap))-gap)/2;
    ymin = -ymax;

    /* hit a slab or a gap ? */

    if (z>zmin && z<zmax && y>ymin && y<ymax) { /* Intersect the crystal? */
      double tilth,tiltv;         /* used to calculate tilt angle of slab */
      double ratio, Q_order, k, kux,kuy,kuz;
      double kix,kiy,kiz;
      int    do_transmit = 0;

      col = ceil ( (z-zmin)/(SlabWidth +gap));  /* which slab hit ? */
      row = ceil ( (y-ymin)/(SlabHeight+gap));
      if (RH != 0) tilth = tiltH ? tiltH[(int)col] :  asin((col-(NH+1)/2)*(SlabWidth+gap)/RH);
      else tilth=0;
      if (RV != 0) tiltv = tiltV ? tiltV[(int)row] : -asin((row-(NV+1)/2)*(SlabHeight+gap)/RV);
      else tiltv=0;

      /* restore neutron in order to transform to slab coordinates */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);

      /* rotate with tilth (around Y) and tiltv (around Z), center on plate */
      double center_z=zmin+(col-0.5)*(SlabWidth+gap) -gap/2;
      double center_y=ymin+(row-0.5)*(SlabHeight+gap)-gap/2;
      Rotation T;
      rot_set_rotation(T, 0, tilth,    tiltv);
      /* now make the coordinate system change */
      mccoordschange_polarisation(T, &vx, &vy, &vz);
      y-=center_y;
      z-=center_z;
      coords_get(rot_apply(T,coords_set(x,y,z)),&x,&y,&z);

      /* this is where polaisation should be handled, plus further down */
      /* mccoordschange_polarisation(t, &sx, &sy, &sz); */

      /* now propagate to slab plane */
      PROP_X0;

      if (fabs(z) <= SlabWidth/2 && fabs(y) <= SlabHeight/2) { /* not in gap ? */
        kix = V2K*vx;             /* Initial wave vector */
        kiy = V2K*vy;
        kiz = V2K*vz;
        /* Get reflection order and corresponding nominal scattering vector q0
          of correct length and direction. Only the order with the closest
          scattering vector is considered */
        ratio = -2*kix/mono_Q;
        Q_order = floor(ratio + .5);
        if(Q_order == 0.0) Q_order = ratio < 0 ? -1 : 1;
        /* Order will be negative when the neutron enters from the back, in
          which case the direction of Q0 is flipped. */
        if(Q_order < 0) Q_order = -Q_order;
        /* Make sure the order is small enough to allow Bragg scattering at the
          given neutron wavelength */
        k = sqrt(kix*kix + kiy*kiy + kiz*kiz);
        kux = kix/k;              /* Unit vector along ki */
        kuy = kiy/k;
        kuz = kiz/k;
        if(Q_order > 2*k/mono_Q) Q_order--;
        if((!order && Q_order > 0) || (Q_order == fabs(order) && order)) {           /* Bragg scattering possible? */
          double q0, q0x, theta, delta, p_reflect, my_r0;

          q0 = Q_order*mono_Q;
          q0x = ratio < 0 ? -q0 : q0;
          theta = asin(q0/(2*k)); /* Actual bragg angle */
          /* Make MC choice: reflect or transmit? */
          delta = asin(fabs(kux)) - theta;

          if (rTable.data != NULL)
          {
            my_r0 = r0*Table_Value(rTable, k, 1); /* 2nd column */
          }
          else my_r0 = r0;
          if (my_r0 > 1)
          {
            if (my_r0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered reflectivity from %f to 1 (k=%f)\n", my_r0, k);
            my_r0=0.999;
          }
          if (my_r0 < 0)
          {
            if (verbose) fprintf(stdout, "Warning: Monochromator_curved : raised reflectivity from %f to 0 (k=%f)\n", my_r0, k);
            my_r0=0;
          }

          p_reflect = fabs(my_r0)*exp(-kiz*kiz/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_y*mos_rms_y))*
                        exp(-kiy*kiy/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_z*mos_rms_z));

          if(rand01() <= p_reflect) { /* Reflect */
            double bx,by,bz,ax,ay,az,phi;
            double cos_2theta,k_sin_2theta,cos_phi,sin_phi,q_x,q_y,q_z;
            double total,c1x,c1y,c1z,w,mos_sample;
            int i=0;

            cos_2theta   = cos(2*theta);
            k_sin_2theta = k*sin(2*theta);
            /* Get unit normal to plane containing ki and most probable kf */
            vec_prod(bx, by, bz, kix, kiy, kiz, q0x, 0, 0);
            NORM(bx,by,bz);
            bx *= k_sin_2theta;
            by *= k_sin_2theta;
            bz *= k_sin_2theta;
            /* Get unit vector normal to ki and b */
            vec_prod(ax, ay, az, bx, by, bz, kux, kuy, kuz);
            /* Compute the total scattering probability at this ki */
            total = 0;
            /* Choose width of Gaussian distribution to sample the angle
            * phi on the Debye-Scherrer cone for the scattered neutron.
            * The radius of the Debye-Scherrer cone is smaller by a
            * factor 1/cos(theta) than the radius of the (partial) sphere
            * describing the possible orientations of Q due to mosaicity, so we
            * start with a width 1/cos(theta) greater than the largest of
            * the two mosaics. */
            mos_sample = mos_rms_max/cos(theta);
            c1x = kix*(cos_2theta-1);
            c1y = kiy*(cos_2theta-1);
            c1z = kiz*(cos_2theta-1);
            /* Loop, repeatedly reducing the sample width until it is small
            * enough to avoid sampling scattering directions with
            * ridiculously low scattering probability.
            * Use a cut-off at 5 times the gauss width for considering
            * scattering probability as well as for integration limits
            * when integrating the sampled distribution below. */
            for(i=0; i<100; i++) {
              w = 5*mos_sample;
              cos_phi = cos(w);
              sin_phi = sin(w);
              q_x =  c1x + cos_phi*ax + sin_phi*bx;
              q_y = (c1y + cos_phi*ay + sin_phi*by)/mos_rms_z;
              q_z = (c1z + cos_phi*az + sin_phi*bz)/mos_rms_y;
              /* Stop when we get near a factor of 25=5^2. */
              if(q_z*q_z + q_y*q_y < (25/(2.0/3.0))*(q_x*q_x))
                break;
              mos_sample *= (2.0/3.0);
            }
            /* Now integrate the chosen sampling distribution, using a
            * cut-off at five times sigma. */
            for(i = 0; i < (sizeof(Gauss_X)/sizeof(double)); i++)
            {
              phi = w*Gauss_X[i];
              cos_phi = cos(phi);
              sin_phi = sin(phi);
              q_x = c1x + cos_phi*ax + sin_phi*bx;
              q_y = c1y + cos_phi*ay + sin_phi*by;
              q_z = c1z + cos_phi*az + sin_phi*bz;
              p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                          GAUSS((q_y/q_x),0,mos_rms_z);
              total += Gauss_W[i]*p_reflect;
            }
            total *= w;
            /* Choose point on Debye-Scherrer cone. Sample from a Gaussian of
             * width 1/cos(theta) greater than the mosaic and correct for any
             * error by adjusting the neutron weight later. */
            phi = mos_sample*randnorm();
            /* Compute final wave vector kf and scattering vector q = ki - kf */
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = c1y + cos_phi*ay + sin_phi*by;
            q_z = c1z + cos_phi*az + sin_phi*bz;
            p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                        GAUSS((q_y/q_x),0,mos_rms_z);

            vx = K2V*(kix+q_x);
            vy = K2V*(kiy+q_y);
            vz = K2V*(kiz+q_z);
            p_reflect /= total*GAUSS(phi,0,mos_sample);
            if (p_reflect <= 0) ABSORB;
            if (p_reflect > 1)  p_reflect = 1;
            p *= p_reflect;

          } /* End MC choice to reflect or transmit neutron (if tmp<p_reflect) */
          else do_transmit = 1;
            /* else transmit neutron */
        } /* End bragg scattering possible (if order) */
        else do_transmit=1;
        if (do_transmit)
        {
          double my_t0;
          if (tTable.data != NULL)
          {
            my_t0 = t0*Table_Value(tTable, k, 1); /* 2nd column */
          }
          else my_t0 = t0;
          /* do not SCATTER, else GROUP does not work */
          if (my_t0 > 1)
          {
            if (my_t0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered transmission from %f to 1 (k=%f)\n", my_t0, k);
            my_t0=0.999;
          }
          if (my_t0 > 0) p*= my_t0;
          else ABSORB;
        }
      } /* end if not in gap */
      /* rotate back in component frame */
      Rotation TT;
      rot_transpose(T, TT);
      /* now make the coordinate system change */
      mccoordschange_polarisation(TT, &vx, &vy, &vz);
      coords_get(rot_apply(TT,coords_set(x,y,z)),&x,&y,&z);
      y+=center_y;
      z+=center_z;
      /* Visualise scattering point in proper, component frame 
	 - but only if the neutron is reflected, that is none of:
	 * transmitted
	 * falling outside the slab material */	
      if(!do_transmit) SCATTER;

      /* mccoordschange_polarisation(tt, &sx, &sy, &sz); */
    } /* End intersect the crystal (if z) */
    else {
      /* restore neutron state when no interaction */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  } /* End neutron moving towards crystal (if vx)*/
}
#line 17539 "./FullInstrument_v4.7_frontend.c"
}   /* End of MonoFLat=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMonoFLat:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(15,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Mono [16] */
  mccoordschange(mcposrMono, mcrotrMono,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Mono (without coords transformations) */
  mcJumpTrace_Mono:
  SIG_MESSAGE("Mono (Trace)");
  mcDEBUG_COMP("Mono")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompMono
  STORE_NEUTRON(16,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[16]++;
  mcPCounter[16] += p;
  mcP2Counter[16] += p*p;
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 16
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define SlabWidth mccMono_SlabWidth
#define SlabHeight mccMono_SlabHeight
#define rTable mccMono_rTable
#define tTable mccMono_tTable
#define row mccMono_row
#define col mccMono_col
#define tiltH mccMono_tiltH
#define tiltV mccMono_tiltV
{   /* Declarations of Mono=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono_reflect;
char* transmit = mccMono_transmit;
MCNUM zwidth = mccMono_zwidth;
MCNUM yheight = mccMono_yheight;
MCNUM gap = mccMono_gap;
MCNUM NH = mccMono_NH;
MCNUM NV = mccMono_NV;
MCNUM mosaich = mccMono_mosaich;
MCNUM mosaicv = mccMono_mosaicv;
MCNUM r0 = mccMono_r0;
MCNUM t0 = mccMono_t0;
MCNUM Q = mccMono_Q;
MCNUM RV = mccMono_RV;
MCNUM RH = mccMono_RH;
MCNUM DM = mccMono_DM;
MCNUM mosaic = mccMono_mosaic;
MCNUM width = mccMono_width;
MCNUM height = mccMono_height;
MCNUM verbose = mccMono_verbose;
MCNUM order = mccMono_order;
/* 'Mono=Monochromator_curved()' component instance has conditional execution */
if (( mcipold_mono == 0 ))

#line 230 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  double dt;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0 && r0)
  {                             /* Moving towards crystal? */
    double zmin,zmax, ymin,ymax;

    /* Propagate to crystal plane */
    PROP_DT(dt);    /* now in the vertical plane of monochromator */

    zmax = ((NH*(SlabWidth+gap))-gap)/2;
    zmin = -zmax;
    ymax = ((NV*(SlabHeight+gap))-gap)/2;
    ymin = -ymax;

    /* hit a slab or a gap ? */

    if (z>zmin && z<zmax && y>ymin && y<ymax) { /* Intersect the crystal? */
      double tilth,tiltv;         /* used to calculate tilt angle of slab */
      double ratio, Q_order, k, kux,kuy,kuz;
      double kix,kiy,kiz;
      int    do_transmit = 0;

      col = ceil ( (z-zmin)/(SlabWidth +gap));  /* which slab hit ? */
      row = ceil ( (y-ymin)/(SlabHeight+gap));
      if (RH != 0) tilth = tiltH ? tiltH[(int)col] :  asin((col-(NH+1)/2)*(SlabWidth+gap)/RH);
      else tilth=0;
      if (RV != 0) tiltv = tiltV ? tiltV[(int)row] : -asin((row-(NV+1)/2)*(SlabHeight+gap)/RV);
      else tiltv=0;

      /* restore neutron in order to transform to slab coordinates */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);

      /* rotate with tilth (around Y) and tiltv (around Z), center on plate */
      double center_z=zmin+(col-0.5)*(SlabWidth+gap) -gap/2;
      double center_y=ymin+(row-0.5)*(SlabHeight+gap)-gap/2;
      Rotation T;
      rot_set_rotation(T, 0, tilth,    tiltv);
      /* now make the coordinate system change */
      mccoordschange_polarisation(T, &vx, &vy, &vz);
      y-=center_y;
      z-=center_z;
      coords_get(rot_apply(T,coords_set(x,y,z)),&x,&y,&z);

      /* this is where polaisation should be handled, plus further down */
      /* mccoordschange_polarisation(t, &sx, &sy, &sz); */

      /* now propagate to slab plane */
      PROP_X0;

      if (fabs(z) <= SlabWidth/2 && fabs(y) <= SlabHeight/2) { /* not in gap ? */
        kix = V2K*vx;             /* Initial wave vector */
        kiy = V2K*vy;
        kiz = V2K*vz;
        /* Get reflection order and corresponding nominal scattering vector q0
          of correct length and direction. Only the order with the closest
          scattering vector is considered */
        ratio = -2*kix/mono_Q;
        Q_order = floor(ratio + .5);
        if(Q_order == 0.0) Q_order = ratio < 0 ? -1 : 1;
        /* Order will be negative when the neutron enters from the back, in
          which case the direction of Q0 is flipped. */
        if(Q_order < 0) Q_order = -Q_order;
        /* Make sure the order is small enough to allow Bragg scattering at the
          given neutron wavelength */
        k = sqrt(kix*kix + kiy*kiy + kiz*kiz);
        kux = kix/k;              /* Unit vector along ki */
        kuy = kiy/k;
        kuz = kiz/k;
        if(Q_order > 2*k/mono_Q) Q_order--;
        if((!order && Q_order > 0) || (Q_order == fabs(order) && order)) {           /* Bragg scattering possible? */
          double q0, q0x, theta, delta, p_reflect, my_r0;

          q0 = Q_order*mono_Q;
          q0x = ratio < 0 ? -q0 : q0;
          theta = asin(q0/(2*k)); /* Actual bragg angle */
          /* Make MC choice: reflect or transmit? */
          delta = asin(fabs(kux)) - theta;

          if (rTable.data != NULL)
          {
            my_r0 = r0*Table_Value(rTable, k, 1); /* 2nd column */
          }
          else my_r0 = r0;
          if (my_r0 > 1)
          {
            if (my_r0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered reflectivity from %f to 1 (k=%f)\n", my_r0, k);
            my_r0=0.999;
          }
          if (my_r0 < 0)
          {
            if (verbose) fprintf(stdout, "Warning: Monochromator_curved : raised reflectivity from %f to 0 (k=%f)\n", my_r0, k);
            my_r0=0;
          }

          p_reflect = fabs(my_r0)*exp(-kiz*kiz/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_y*mos_rms_y))*
                        exp(-kiy*kiy/(kiy*kiy + kiz*kiz)*(delta*delta)/
                            (2*mos_rms_z*mos_rms_z));

          if(rand01() <= p_reflect) { /* Reflect */
            double bx,by,bz,ax,ay,az,phi;
            double cos_2theta,k_sin_2theta,cos_phi,sin_phi,q_x,q_y,q_z;
            double total,c1x,c1y,c1z,w,mos_sample;
            int i=0;

            cos_2theta   = cos(2*theta);
            k_sin_2theta = k*sin(2*theta);
            /* Get unit normal to plane containing ki and most probable kf */
            vec_prod(bx, by, bz, kix, kiy, kiz, q0x, 0, 0);
            NORM(bx,by,bz);
            bx *= k_sin_2theta;
            by *= k_sin_2theta;
            bz *= k_sin_2theta;
            /* Get unit vector normal to ki and b */
            vec_prod(ax, ay, az, bx, by, bz, kux, kuy, kuz);
            /* Compute the total scattering probability at this ki */
            total = 0;
            /* Choose width of Gaussian distribution to sample the angle
            * phi on the Debye-Scherrer cone for the scattered neutron.
            * The radius of the Debye-Scherrer cone is smaller by a
            * factor 1/cos(theta) than the radius of the (partial) sphere
            * describing the possible orientations of Q due to mosaicity, so we
            * start with a width 1/cos(theta) greater than the largest of
            * the two mosaics. */
            mos_sample = mos_rms_max/cos(theta);
            c1x = kix*(cos_2theta-1);
            c1y = kiy*(cos_2theta-1);
            c1z = kiz*(cos_2theta-1);
            /* Loop, repeatedly reducing the sample width until it is small
            * enough to avoid sampling scattering directions with
            * ridiculously low scattering probability.
            * Use a cut-off at 5 times the gauss width for considering
            * scattering probability as well as for integration limits
            * when integrating the sampled distribution below. */
            for(i=0; i<100; i++) {
              w = 5*mos_sample;
              cos_phi = cos(w);
              sin_phi = sin(w);
              q_x =  c1x + cos_phi*ax + sin_phi*bx;
              q_y = (c1y + cos_phi*ay + sin_phi*by)/mos_rms_z;
              q_z = (c1z + cos_phi*az + sin_phi*bz)/mos_rms_y;
              /* Stop when we get near a factor of 25=5^2. */
              if(q_z*q_z + q_y*q_y < (25/(2.0/3.0))*(q_x*q_x))
                break;
              mos_sample *= (2.0/3.0);
            }
            /* Now integrate the chosen sampling distribution, using a
            * cut-off at five times sigma. */
            for(i = 0; i < (sizeof(Gauss_X)/sizeof(double)); i++)
            {
              phi = w*Gauss_X[i];
              cos_phi = cos(phi);
              sin_phi = sin(phi);
              q_x = c1x + cos_phi*ax + sin_phi*bx;
              q_y = c1y + cos_phi*ay + sin_phi*by;
              q_z = c1z + cos_phi*az + sin_phi*bz;
              p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                          GAUSS((q_y/q_x),0,mos_rms_z);
              total += Gauss_W[i]*p_reflect;
            }
            total *= w;
            /* Choose point on Debye-Scherrer cone. Sample from a Gaussian of
             * width 1/cos(theta) greater than the mosaic and correct for any
             * error by adjusting the neutron weight later. */
            phi = mos_sample*randnorm();
            /* Compute final wave vector kf and scattering vector q = ki - kf */
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = c1y + cos_phi*ay + sin_phi*by;
            q_z = c1z + cos_phi*az + sin_phi*bz;
            p_reflect = GAUSS((q_z/q_x),0,mos_rms_y)*
                        GAUSS((q_y/q_x),0,mos_rms_z);

            vx = K2V*(kix+q_x);
            vy = K2V*(kiy+q_y);
            vz = K2V*(kiz+q_z);
            p_reflect /= total*GAUSS(phi,0,mos_sample);
            if (p_reflect <= 0) ABSORB;
            if (p_reflect > 1)  p_reflect = 1;
            p *= p_reflect;

          } /* End MC choice to reflect or transmit neutron (if tmp<p_reflect) */
          else do_transmit = 1;
            /* else transmit neutron */
        } /* End bragg scattering possible (if order) */
        else do_transmit=1;
        if (do_transmit)
        {
          double my_t0;
          if (tTable.data != NULL)
          {
            my_t0 = t0*Table_Value(tTable, k, 1); /* 2nd column */
          }
          else my_t0 = t0;
          /* do not SCATTER, else GROUP does not work */
          if (my_t0 > 1)
          {
            if (my_t0 > 1.01 && verbose) fprintf(stdout, "Warning: Monochromator_curved : lowered transmission from %f to 1 (k=%f)\n", my_t0, k);
            my_t0=0.999;
          }
          if (my_t0 > 0) p*= my_t0;
          else ABSORB;
        }
      } /* end if not in gap */
      /* rotate back in component frame */
      Rotation TT;
      rot_transpose(T, TT);
      /* now make the coordinate system change */
      mccoordschange_polarisation(TT, &vx, &vy, &vz);
      coords_get(rot_apply(TT,coords_set(x,y,z)),&x,&y,&z);
      y+=center_y;
      z+=center_z;
      /* Visualise scattering point in proper, component frame 
	 - but only if the neutron is reflected, that is none of:
	 * transmitted
	 * falling outside the slab material */	
      if(!do_transmit) SCATTER;

      /* mccoordschange_polarisation(tt, &sx, &sy, &sz); */
    } /* End intersect the crystal (if z) */
    else {
      /* restore neutron state when no interaction */
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
  } /* End neutron moving towards crystal (if vx)*/
}
#line 17920 "./FullInstrument_v4.7_frontend.c"
}   /* End of Mono=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompMono:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(16,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component A2 [17] */
  mccoordschange(mcposrA2, mcrotrA2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component A2 (without coords transformations) */
  mcJumpTrace_A2:
  SIG_MESSAGE("A2 (Trace)");
  mcDEBUG_COMP("A2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompA2
  STORE_NEUTRON(17,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[17]++;
  mcPCounter[17] += p;
  mcP2Counter[17] += p*p;
#define mccompcurname  A2
#define mccompcurtype  Arm
#define mccompcurindex 17
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompA2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(17,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component slit [18] */
  mccoordschange(mcposrslit, mcrotrslit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit (without coords transformations) */
  mcJumpTrace_slit:
  SIG_MESSAGE("slit (Trace)");
  mcDEBUG_COMP("slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompslit
  STORE_NEUTRON(18,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[18]++;
  mcPCounter[18] += p;
  mcP2Counter[18] += p*p;
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of slit=Slit() SETTING parameters. */
MCNUM xmin = mccslit_xmin;
MCNUM xmax = mccslit_xmax;
MCNUM ymin = mccslit_ymin;
MCNUM ymax = mccslit_ymax;
MCNUM radius = mccslit_radius;
MCNUM xwidth = mccslit_xwidth;
MCNUM yheight = mccslit_yheight;
#line 71 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 18159 "./FullInstrument_v4.7_frontend.c"
}   /* End of slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(18,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component SamplePos [19] */
  mccoordschange(mcposrSamplePos, mcrotrSamplePos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component SamplePos (without coords transformations) */
  mcJumpTrace_SamplePos:
  SIG_MESSAGE("SamplePos (Trace)");
  mcDEBUG_COMP("SamplePos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSamplePos
  STORE_NEUTRON(19,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[19]++;
  mcPCounter[19] += p;
  mcP2Counter[19] += p*p;
#define mccompcurname  SamplePos
#define mccompcurtype  Arm
#define mccompcurindex 19
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSamplePos:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(19,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_sample_pos [20] */
  mccoordschange(mcposrpsd_monitor_sample_pos, mcrotrpsd_monitor_sample_pos,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_sample_pos (without coords transformations) */
  mcJumpTrace_psd_monitor_sample_pos:
  SIG_MESSAGE("psd_monitor_sample_pos (Trace)");
  mcDEBUG_COMP("psd_monitor_sample_pos")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_sample_pos
  STORE_NEUTRON(20,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[20]++;
  mcPCounter[20] += p;
  mcP2Counter[20] += p*p;
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
{   /* Declarations of psd_monitor_sample_pos=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_sample_pos_nx;
int ny = mccpsd_monitor_sample_pos_ny;
char* filename = mccpsd_monitor_sample_pos_filename;
MCNUM xmin = mccpsd_monitor_sample_pos_xmin;
MCNUM xmax = mccpsd_monitor_sample_pos_xmax;
MCNUM ymin = mccpsd_monitor_sample_pos_ymin;
MCNUM ymax = mccpsd_monitor_sample_pos_ymax;
MCNUM xwidth = mccpsd_monitor_sample_pos_xwidth;
MCNUM yheight = mccpsd_monitor_sample_pos_yheight;
MCNUM restore_neutron = mccpsd_monitor_sample_pos_restore_neutron;
int nowritefile = mccpsd_monitor_sample_pos_nowritefile;
#line 95 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 18397 "./FullInstrument_v4.7_frontend.c"
}   /* End of psd_monitor_sample_pos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_sample_pos:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(20,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component e_monitor_sample_ps [21] */
  mccoordschange(mcposre_monitor_sample_ps, mcrotre_monitor_sample_ps,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component e_monitor_sample_ps (without coords transformations) */
  mcJumpTrace_e_monitor_sample_ps:
  SIG_MESSAGE("e_monitor_sample_ps (Trace)");
  mcDEBUG_COMP("e_monitor_sample_ps")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompe_monitor_sample_ps
  STORE_NEUTRON(21,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[21]++;
  mcPCounter[21] += p;
  mcP2Counter[21] += p*p;
#define mccompcurname  e_monitor_sample_ps
#define mccompcurtype  E_monitor
#define mccompcurindex 21
#define nE mcce_monitor_sample_ps_nE
#define E_N mcce_monitor_sample_ps_E_N
#define E_p mcce_monitor_sample_ps_E_p
#define E_p2 mcce_monitor_sample_ps_E_p2
#define S_p mcce_monitor_sample_ps_S_p
#define S_pE mcce_monitor_sample_ps_S_pE
#define S_pE2 mcce_monitor_sample_ps_S_pE2
{   /* Declarations of e_monitor_sample_ps=E_monitor() SETTING parameters. */
char* filename = mcce_monitor_sample_ps_filename;
MCNUM xmin = mcce_monitor_sample_ps_xmin;
MCNUM xmax = mcce_monitor_sample_ps_xmax;
MCNUM ymin = mcce_monitor_sample_ps_ymin;
MCNUM ymax = mcce_monitor_sample_ps_ymax;
MCNUM xwidth = mcce_monitor_sample_ps_xwidth;
MCNUM yheight = mcce_monitor_sample_ps_yheight;
MCNUM Emin = mcce_monitor_sample_ps_Emin;
MCNUM Emax = mcce_monitor_sample_ps_Emax;
MCNUM restore_neutron = mcce_monitor_sample_ps_restore_neutron;
int nowritefile = mcce_monitor_sample_ps_nowritefile;
#line 89 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
    int i;
    double E;

    PROP_Z0;
    if (x>xmin && x<xmax && y>ymin && y<ymax)
    {
      E = VS2E*(vx*vx + vy*vy + vz*vz);

      S_p += p;
      S_pE += p*E;
      S_pE2 += p*E*E;

      i = floor((E-Emin)*nE/(Emax-Emin));
      if(i >= 0 && i < nE)
      {
        E_N[i]++;
        E_p[i] += p;
        E_p2[i] += p*p;
        SCATTER;
      }
    }
    if (restore_neutron) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
    }
}
#line 18551 "./FullInstrument_v4.7_frontend.c"
}   /* End of e_monitor_sample_ps=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompe_monitor_sample_ps:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(21,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'Origin'. */
  SIG_MESSAGE("Origin (Save)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 115 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 18667 "./FullInstrument_v4.7_frontend.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'SlitMonitor'. */
  SIG_MESSAGE("SlitMonitor (Save)");
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
{   /* Declarations of SlitMonitor=PSD_monitor() SETTING parameters. */
int nx = mccSlitMonitor_nx;
int ny = mccSlitMonitor_ny;
char* filename = mccSlitMonitor_filename;
MCNUM xmin = mccSlitMonitor_xmin;
MCNUM xmax = mccSlitMonitor_xmax;
MCNUM ymin = mccSlitMonitor_ymin;
MCNUM ymax = mccSlitMonitor_ymax;
MCNUM xwidth = mccSlitMonitor_xwidth;
MCNUM yheight = mccSlitMonitor_yheight;
MCNUM restore_neutron = mccSlitMonitor_restore_neutron;
int nowritefile = mccSlitMonitor_nowritefile;
#line 111 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
if (!nowritefile) {
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
}
#line 18710 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitMonitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_sample_pos'. */
  SIG_MESSAGE("psd_monitor_sample_pos (Save)");
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
{   /* Declarations of psd_monitor_sample_pos=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_sample_pos_nx;
int ny = mccpsd_monitor_sample_pos_ny;
char* filename = mccpsd_monitor_sample_pos_filename;
MCNUM xmin = mccpsd_monitor_sample_pos_xmin;
MCNUM xmax = mccpsd_monitor_sample_pos_xmax;
MCNUM ymin = mccpsd_monitor_sample_pos_ymin;
MCNUM ymax = mccpsd_monitor_sample_pos_ymax;
MCNUM xwidth = mccpsd_monitor_sample_pos_xwidth;
MCNUM yheight = mccpsd_monitor_sample_pos_yheight;
MCNUM restore_neutron = mccpsd_monitor_sample_pos_restore_neutron;
int nowritefile = mccpsd_monitor_sample_pos_nowritefile;
#line 111 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
if (!nowritefile) {
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
}
#line 18752 "./FullInstrument_v4.7_frontend.c"
}   /* End of psd_monitor_sample_pos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'e_monitor_sample_ps'. */
  SIG_MESSAGE("e_monitor_sample_ps (Save)");
#define mccompcurname  e_monitor_sample_ps
#define mccompcurtype  E_monitor
#define mccompcurindex 21
#define nE mcce_monitor_sample_ps_nE
#define E_N mcce_monitor_sample_ps_E_N
#define E_p mcce_monitor_sample_ps_E_p
#define E_p2 mcce_monitor_sample_ps_E_p2
#define S_p mcce_monitor_sample_ps_S_p
#define S_pE mcce_monitor_sample_ps_S_pE
#define S_pE2 mcce_monitor_sample_ps_S_pE2
{   /* Declarations of e_monitor_sample_ps=E_monitor() SETTING parameters. */
char* filename = mcce_monitor_sample_ps_filename;
MCNUM xmin = mcce_monitor_sample_ps_xmin;
MCNUM xmax = mcce_monitor_sample_ps_xmax;
MCNUM ymin = mcce_monitor_sample_ps_ymin;
MCNUM ymax = mcce_monitor_sample_ps_ymax;
MCNUM xwidth = mcce_monitor_sample_ps_xwidth;
MCNUM yheight = mcce_monitor_sample_ps_yheight;
MCNUM Emin = mcce_monitor_sample_ps_Emin;
MCNUM Emax = mcce_monitor_sample_ps_Emax;
MCNUM restore_neutron = mcce_monitor_sample_ps_restore_neutron;
int nowritefile = mcce_monitor_sample_ps_nowritefile;
#line 117 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
    if (!nowritefile) {
    DETECTOR_OUT_1D(
        "Energy monitor",
        "Energy [meV]",
        "Intensity",
        "E", Emin, Emax, nE,
        &E_N[0],&E_p[0],&E_p2[0],
        filename);
    if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
     S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
    }
}
#line 18799 "./FullInstrument_v4.7_frontend.c"
}   /* End of e_monitor_sample_ps=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'Origin'. */
  SIG_MESSAGE("Origin (Finally)");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 133 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 18846 "./FullInstrument_v4.7_frontend.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
  /* User FINALLY code for component 'source'. */
  SIG_MESSAGE("source (Finally)");
#define mccompcurname  source
#define mccompcurtype  Source_gen4
#define mccompcurindex 2
#define flux_file mccsource_flux_file
#define xdiv_file mccsource_xdiv_file
#define ydiv_file mccsource_ydiv_file
#define p_in mccsource_p_in
#define lambda0 mccsource_lambda0
#define lambda02 mccsource_lambda02
#define L2P mccsource_L2P
#define lambda0b mccsource_lambda0b
#define lambda02b mccsource_lambda02b
#define L2Pb mccsource_L2Pb
#define lambda0c mccsource_lambda0c
#define lambda02c mccsource_lambda02c
#define L2Pc mccsource_L2Pc
#define pTable mccsource_pTable
#define pTable_x mccsource_pTable_x
#define pTable_y mccsource_pTable_y
{   /* Declarations of source=Source_gen4() SETTING parameters. */
MCNUM radius = mccsource_radius;
MCNUM dist = mccsource_dist;
MCNUM xw = mccsource_xw;
MCNUM yh = mccsource_yh;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM Lambda0 = mccsource_Lambda0;
MCNUM dLambda = mccsource_dLambda;
MCNUM I1 = mccsource_I1;
MCNUM h = mccsource_h;
MCNUM w = mccsource_w;
MCNUM gaussian = mccsource_gaussian;
MCNUM verbose = mccsource_verbose;
MCNUM T1 = mccsource_T1;
MCNUM flux_file_perAA = mccsource_flux_file_perAA;
MCNUM flux_file_log = mccsource_flux_file_log;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM Emin = mccsource_Emin;
MCNUM Emax = mccsource_Emax;
MCNUM T2 = mccsource_T2;
MCNUM I2 = mccsource_I2;
MCNUM T3 = mccsource_T3;
MCNUM I3 = mccsource_I3;
MCNUM length = mccsource_length;
MCNUM phi_init = mccsource_phi_init;
MCNUM theta_init = mccsource_theta_init;
MCNUM HEtailA = mccsource_HEtailA;
MCNUM HEtailL0 = mccsource_HEtailL0;
#line 564 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
{
  Table_Free(&pTable);
  Table_Free(&pTable_x);
  Table_Free(&pTable_y);
}
#line 18915 "./FullInstrument_v4.7_frontend.c"
}   /* End of source=Source_gen4() SETTING parameter declarations. */
#undef pTable_y
#undef pTable_x
#undef pTable
#undef L2Pc
#undef lambda02c
#undef lambda0c
#undef L2Pb
#undef lambda02b
#undef lambda0b
#undef L2P
#undef lambda02
#undef lambda0
#undef p_in
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source=Source_gen4()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] StartOfGuide\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] StartOfGuide=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] elliptical_guide_gravity5\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] elliptical_guide_gravity5=Elliptic_guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] EndOfelement_5\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] EndOfelement_5=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] EndOfelement_4\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] EndOfelement_4=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] elliptical_guide_gravity3\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] elliptical_guide_gravity3=Elliptic_guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] EndOfelement_3\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] EndOfelement_3=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] Curved_guide_2\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] Curved_guide_2=Bender()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] EndOfelement_2\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] EndOfelement_2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] elliptical_guide_gravity1\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] elliptical_guide_gravity1=Elliptic_guide_gravity()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] EndOfelement_1\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] EndOfelement_1=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] SlitGuideEnd\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] SlitGuideEnd=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
  /* User FINALLY code for component 'SlitMonitor'. */
  SIG_MESSAGE("SlitMonitor (Finally)");
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
{   /* Declarations of SlitMonitor=PSD_monitor() SETTING parameters. */
int nx = mccSlitMonitor_nx;
int ny = mccSlitMonitor_ny;
char* filename = mccSlitMonitor_filename;
MCNUM xmin = mccSlitMonitor_xmin;
MCNUM xmax = mccSlitMonitor_xmax;
MCNUM ymin = mccSlitMonitor_ymin;
MCNUM ymax = mccSlitMonitor_ymax;
MCNUM xwidth = mccSlitMonitor_xwidth;
MCNUM yheight = mccSlitMonitor_yheight;
MCNUM restore_neutron = mccSlitMonitor_restore_neutron;
int nowritefile = mccSlitMonitor_nowritefile;
#line 125 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 18987 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitMonitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] SlitMonitor\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] SlitMonitor=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
  /* User FINALLY code for component 'MonoFLat'. */
  SIG_MESSAGE("MonoFLat (Finally)");
#define mccompcurname  MonoFLat
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 15
#define mos_rms_y mccMonoFLat_mos_rms_y
#define mos_rms_z mccMonoFLat_mos_rms_z
#define mos_rms_max mccMonoFLat_mos_rms_max
#define mono_Q mccMonoFLat_mono_Q
#define SlabWidth mccMonoFLat_SlabWidth
#define SlabHeight mccMonoFLat_SlabHeight
#define rTable mccMonoFLat_rTable
#define tTable mccMonoFLat_tTable
#define row mccMonoFLat_row
#define col mccMonoFLat_col
#define tiltH mccMonoFLat_tiltH
#define tiltV mccMonoFLat_tiltV
{   /* Declarations of MonoFLat=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonoFLat_reflect;
char* transmit = mccMonoFLat_transmit;
MCNUM zwidth = mccMonoFLat_zwidth;
MCNUM yheight = mccMonoFLat_yheight;
MCNUM gap = mccMonoFLat_gap;
MCNUM NH = mccMonoFLat_NH;
MCNUM NV = mccMonoFLat_NV;
MCNUM mosaich = mccMonoFLat_mosaich;
MCNUM mosaicv = mccMonoFLat_mosaicv;
MCNUM r0 = mccMonoFLat_r0;
MCNUM t0 = mccMonoFLat_t0;
MCNUM Q = mccMonoFLat_Q;
MCNUM RV = mccMonoFLat_RV;
MCNUM RH = mccMonoFLat_RH;
MCNUM DM = mccMonoFLat_DM;
MCNUM mosaic = mccMonoFLat_mosaic;
MCNUM width = mccMonoFLat_width;
MCNUM height = mccMonoFLat_height;
MCNUM verbose = mccMonoFLat_verbose;
MCNUM order = mccMonoFLat_order;
#line 460 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  Table_Free(&rTable);
  Table_Free(&tTable);
  if (tiltH) free(tiltH);
  if (tiltV) free(tiltV);
}
#line 19043 "./FullInstrument_v4.7_frontend.c"
}   /* End of MonoFLat=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] MonoFLat\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] MonoFLat=Monochromator_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
  /* User FINALLY code for component 'Mono'. */
  SIG_MESSAGE("Mono (Finally)");
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 16
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define SlabWidth mccMono_SlabWidth
#define SlabHeight mccMono_SlabHeight
#define rTable mccMono_rTable
#define tTable mccMono_tTable
#define row mccMono_row
#define col mccMono_col
#define tiltH mccMono_tiltH
#define tiltV mccMono_tiltV
{   /* Declarations of Mono=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono_reflect;
char* transmit = mccMono_transmit;
MCNUM zwidth = mccMono_zwidth;
MCNUM yheight = mccMono_yheight;
MCNUM gap = mccMono_gap;
MCNUM NH = mccMono_NH;
MCNUM NV = mccMono_NV;
MCNUM mosaich = mccMono_mosaich;
MCNUM mosaicv = mccMono_mosaicv;
MCNUM r0 = mccMono_r0;
MCNUM t0 = mccMono_t0;
MCNUM Q = mccMono_Q;
MCNUM RV = mccMono_RV;
MCNUM RH = mccMono_RH;
MCNUM DM = mccMono_DM;
MCNUM mosaic = mccMono_mosaic;
MCNUM width = mccMono_width;
MCNUM height = mccMono_height;
MCNUM verbose = mccMono_verbose;
MCNUM order = mccMono_order;
#line 460 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  Table_Free(&rTable);
  Table_Free(&tTable);
  if (tiltH) free(tiltH);
  if (tiltV) free(tiltV);
}
#line 19108 "./FullInstrument_v4.7_frontend.c"
}   /* End of Mono=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] Mono\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] Mono=Monochromator_curved()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] A2\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] A2=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] slit\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] SamplePos\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] SamplePos=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
  /* User FINALLY code for component 'psd_monitor_sample_pos'. */
  SIG_MESSAGE("psd_monitor_sample_pos (Finally)");
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
{   /* Declarations of psd_monitor_sample_pos=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_sample_pos_nx;
int ny = mccpsd_monitor_sample_pos_ny;
char* filename = mccpsd_monitor_sample_pos_filename;
MCNUM xmin = mccpsd_monitor_sample_pos_xmin;
MCNUM xmax = mccpsd_monitor_sample_pos_xmax;
MCNUM ymin = mccpsd_monitor_sample_pos_ymin;
MCNUM ymax = mccpsd_monitor_sample_pos_ymax;
MCNUM xwidth = mccpsd_monitor_sample_pos_xwidth;
MCNUM yheight = mccpsd_monitor_sample_pos_yheight;
MCNUM restore_neutron = mccpsd_monitor_sample_pos_restore_neutron;
int nowritefile = mccpsd_monitor_sample_pos_nowritefile;
#line 125 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 19160 "./FullInstrument_v4.7_frontend.c"
}   /* End of psd_monitor_sample_pos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] psd_monitor_sample_pos\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] psd_monitor_sample_pos=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] e_monitor_sample_ps\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] e_monitor_sample_ps=E_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
#define cylinder mcdis_cylinder
#define sphere mcdis_sphere
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'Origin'. */
  SIG_MESSAGE("Origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Origin");
#define mccompcurname  Origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccOrigin_IntermediateCnts
#define StartTime mccOrigin_StartTime
#define EndTime mccOrigin_EndTime
#define CurrentTime mccOrigin_CurrentTime
{   /* Declarations of Origin=Progress_bar() SETTING parameters. */
char* profile = mccOrigin_profile;
MCNUM percent = mccOrigin_percent;
MCNUM flag_save = mccOrigin_flag_save;
MCNUM minutes = mccOrigin_minutes;
#line 147 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 19207 "./FullInstrument_v4.7_frontend.c"
}   /* End of Origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source'. */
  SIG_MESSAGE("source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source");
#define mccompcurname  source
#define mccompcurtype  Source_gen4
#define mccompcurindex 2
#define flux_file mccsource_flux_file
#define xdiv_file mccsource_xdiv_file
#define ydiv_file mccsource_ydiv_file
#define p_in mccsource_p_in
#define lambda0 mccsource_lambda0
#define lambda02 mccsource_lambda02
#define L2P mccsource_L2P
#define lambda0b mccsource_lambda0b
#define lambda02b mccsource_lambda02b
#define L2Pb mccsource_L2Pb
#define lambda0c mccsource_lambda0c
#define lambda02c mccsource_lambda02c
#define L2Pc mccsource_L2Pc
#define pTable mccsource_pTable
#define pTable_x mccsource_pTable_x
#define pTable_y mccsource_pTable_y
{   /* Declarations of source=Source_gen4() SETTING parameters. */
MCNUM radius = mccsource_radius;
MCNUM dist = mccsource_dist;
MCNUM xw = mccsource_xw;
MCNUM yh = mccsource_yh;
MCNUM E0 = mccsource_E0;
MCNUM dE = mccsource_dE;
MCNUM Lambda0 = mccsource_Lambda0;
MCNUM dLambda = mccsource_dLambda;
MCNUM I1 = mccsource_I1;
MCNUM h = mccsource_h;
MCNUM w = mccsource_w;
MCNUM gaussian = mccsource_gaussian;
MCNUM verbose = mccsource_verbose;
MCNUM T1 = mccsource_T1;
MCNUM flux_file_perAA = mccsource_flux_file_perAA;
MCNUM flux_file_log = mccsource_flux_file_log;
MCNUM Lmin = mccsource_Lmin;
MCNUM Lmax = mccsource_Lmax;
MCNUM Emin = mccsource_Emin;
MCNUM Emax = mccsource_Emax;
MCNUM T2 = mccsource_T2;
MCNUM I2 = mccsource_I2;
MCNUM T3 = mccsource_T3;
MCNUM I3 = mccsource_I3;
MCNUM length = mccsource_length;
MCNUM phi_init = mccsource_phi_init;
MCNUM theta_init = mccsource_theta_init;
MCNUM HEtailA = mccsource_HEtailA;
MCNUM HEtailL0 = mccsource_HEtailL0;
#line 571 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../contrib/Source_gen4.comp"
{
  double xmin;
  double xmax;
  double ymin;
  double ymax;

  if ((h == 0) || (w == 0))
  {
    
    circle("xy",0,0,0,radius);
    if (gaussian)
      circle("xy",0,0,0,radius/2);
  }
  else
  {
    xmin = -w/2; xmax = w/2;
    ymin = -h/2; ymax = h/2;

    
    multiline(5, (double)xmin, (double)ymin, 0.0,
             (double)xmax, (double)ymin, 0.0,
             (double)xmax, (double)ymax, 0.0,
             (double)xmin, (double)ymax, 0.0,
             (double)xmin, (double)ymin, 0.0);
    if (gaussian)
      circle("xy",0,0,0,sqrt(w*w+h*h)/4);
  }
}
#line 19298 "./FullInstrument_v4.7_frontend.c"
}   /* End of source=Source_gen4() SETTING parameter declarations. */
#undef pTable_y
#undef pTable_x
#undef pTable
#undef L2Pc
#undef lambda02c
#undef lambda0c
#undef L2Pb
#undef lambda02b
#undef lambda0b
#undef L2P
#undef lambda02
#undef lambda0
#undef p_in
#undef ydiv_file
#undef xdiv_file
#undef flux_file
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'StartOfGuide'. */
  SIG_MESSAGE("StartOfGuide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "StartOfGuide");
#define mccompcurname  StartOfGuide
#define mccompcurtype  Arm
#define mccompcurindex 3
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19334 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'elliptical_guide_gravity5'. */
  SIG_MESSAGE("elliptical_guide_gravity5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "elliptical_guide_gravity5");
#define mccompcurname  elliptical_guide_gravity5
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 4
#define mvaluesright mccelliptical_guide_gravity5_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity5_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity5_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity5_mvaluesbottom
#define seglength mccelliptical_guide_gravity5_seglength
#define guideInfo mccelliptical_guide_gravity5_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity5_latestParticleCollision
#define Gx mccelliptical_guide_gravity5_Gx
#define Gy mccelliptical_guide_gravity5_Gy
#define Gz mccelliptical_guide_gravity5_Gz
#define Gx0 mccelliptical_guide_gravity5_Gx0
#define Gy0 mccelliptical_guide_gravity5_Gy0
#define Gz0 mccelliptical_guide_gravity5_Gz0
#define Circ mccelliptical_guide_gravity5_Circ
#define dynamicalSegLength mccelliptical_guide_gravity5_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity5=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity5_l;
MCNUM xwidth = mccelliptical_guide_gravity5_xwidth;
MCNUM yheight = mccelliptical_guide_gravity5_yheight;
MCNUM linxw = mccelliptical_guide_gravity5_linxw;
MCNUM loutxw = mccelliptical_guide_gravity5_loutxw;
MCNUM linyh = mccelliptical_guide_gravity5_linyh;
MCNUM loutyh = mccelliptical_guide_gravity5_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity5_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity5_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity5_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity5_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity5_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity5_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity5_dimensionsAt;
char* option = mccelliptical_guide_gravity5_option;
MCNUM R0 = mccelliptical_guide_gravity5_R0;
MCNUM Qc = mccelliptical_guide_gravity5_Qc;
MCNUM alpha = mccelliptical_guide_gravity5_alpha;
MCNUM m = mccelliptical_guide_gravity5_m;
MCNUM W = mccelliptical_guide_gravity5_W;
MCNUM alpharight = mccelliptical_guide_gravity5_alpharight;
MCNUM mright = mccelliptical_guide_gravity5_mright;
MCNUM alphaleft = mccelliptical_guide_gravity5_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity5_mleft;
MCNUM alphatop = mccelliptical_guide_gravity5_alphatop;
MCNUM mtop = mccelliptical_guide_gravity5_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity5_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity5_mbottom;
char* verbose = mccelliptical_guide_gravity5_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity5_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity5_curvature;
#line 1496 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{


	// Calculate the points need to draw approximation of the ellipses
	// defining the guide 

	// the number of lines used to draw one side of the guide
	int ApproximationMirrors = 500; 

	// The start of the guide
	double zvalue=0;

	// The the different in z between point used to draw the ellipse
	double zdelta = guideInfo.Length/(1.0*ApproximationMirrors);

	// The vector used to store the points defining the lines
	double xplus[ApproximationMirrors+1];
	double xminus[ApproximationMirrors+1];
	double yplus[ApproximationMirrors+1];
	double yminus[ApproximationMirrors+1];

	// Temperary values for the loop
	double tempx;
	double tempy;

	/*
		Calculate the second coordinates to the points on the ellipse with z_i 
		as the first coordinate. We transform the point to the coordinate system 
		there the ellipse is a unit circle. And use the defination of this circle 
		to find second coordinate (x^2+z^2 = 1)
	*/

	/////////////////////////////////////////////////////////
	double Length;
	double entranceHorizontalWidth;
	double entranceVerticalWidth;	
	
	// ellipses infomation
	double ellipseMajorAxis[4],		ellipseMinorAxis[4];
	double ellipseMajorOffset[4],	ellipseMinorOffset[4];

	enum Side {RightSide,TopSide,LeftSide,BottomSide,None};
	/////////////////////////////////////////////////////////

	int i = 0;
	double tempz = 0;
	for(i=0;i<ApproximationMirrors+1;i++){
		
		tempx = sqrt(
			guideInfo.ellipseMinorAxis[RightSide]
			*guideInfo.ellipseMinorAxis[RightSide]
			-(	guideInfo.ellipseMinorAxis[RightSide]
				*guideInfo.ellipseMinorAxis[RightSide] )
			/(	guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
		);

		xplus[i] =  tempx + guideInfo.ellipseMinorOffset[RightSide];
		xminus[i]= -tempx + guideInfo.ellipseMinorOffset[RightSide];
		
		tempy = sqrt(
			guideInfo.ellipseMinorAxis[TopSide]
			*guideInfo.ellipseMinorAxis[TopSide]
			-(	guideInfo.ellipseMinorAxis[TopSide]
				*guideInfo.ellipseMinorAxis[TopSide] )
			/(	guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
		);
		
		yplus[i] =  tempy + guideInfo.ellipseMinorOffset[TopSide];
		yminus[i]= -tempy + guideInfo.ellipseMinorOffset[TopSide];
	}

	///// Draw lines

	// Drawing lines orthogonal with the z direction.
	// at both ends of the guide and at the boardest place at the guide

	// These may not give correct result if one of the ends are closed

	int j=0;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 	  , yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta,xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j],0 		, zvalue+j*zdelta);
			
	j=ApproximationMirrors;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j],zvalue+j*zdelta,xminus[j] , yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	    , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		 , zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j], 0 		 , zvalue+j*zdelta);

			  
			  
	// find boardest place on the guide and draw a band around the guide
	int m0;
	double boardestPlace = 0;
	int boardestPlaceNumber = 0;
	for(m0=0; m0<ApproximationMirrors; m0++){
		if( boardestPlace <= fabs(yplus[m0]) ){
			boardestPlace = fabs(yplus[m0]);
			boardestPlaceNumber = m0;
		}
	}
	j = boardestPlaceNumber;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta, xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta, xplus[j], 0 		 , zvalue+j*zdelta);


			  
	// Drawing lines parallel with the z direction
			  
	int k=0;
	for(k=0; k < ApproximationMirrors; k++){
		
		line( xplus[k], yplus[k], zvalue+k*zdelta,xplus[k+1], yplus[k+1], zvalue+(k+1)*zdelta);
		line( xminus[k],yplus[k], zvalue+k*zdelta,xminus[k+1],yplus[k+1], zvalue+(k+1)*zdelta);
		
		line( xplus[k], yminus[k],zvalue+k*zdelta,xplus[k+1], yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],yminus[k],zvalue+k*zdelta,xminus[k+1],yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],0 		, zvalue+k*zdelta, xminus[k+1],0 	, zvalue+(k+1)*zdelta);
		line( xplus[k], 0 		, zvalue+k*zdelta, xplus[k+1], 0 	, zvalue+(k+1)*zdelta);
		
		line( 0 	, yminus[k],  zvalue+k*zdelta, 0 	  	,yminus[k+1],zvalue+(k+1)*zdelta);
		line( 0 	,yplus[k]  , zvalue+k*zdelta , 0 		,yplus[k] , zvalue+(k+1)*zdelta);
		
		}

		if(guideInfo.EnclosingBoxOn){
			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );

			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );
			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
		}
}
#line 19582 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity5=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EndOfelement_5'. */
  SIG_MESSAGE("EndOfelement_5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EndOfelement_5");
#define mccompcurname  EndOfelement_5
#define mccompcurtype  Arm
#define mccompcurindex 5
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19617 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EndOfelement_4'. */
  SIG_MESSAGE("EndOfelement_4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EndOfelement_4");
#define mccompcurname  EndOfelement_4
#define mccompcurtype  Arm
#define mccompcurindex 6
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19636 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'elliptical_guide_gravity3'. */
  SIG_MESSAGE("elliptical_guide_gravity3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "elliptical_guide_gravity3");
#define mccompcurname  elliptical_guide_gravity3
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 7
#define mvaluesright mccelliptical_guide_gravity3_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity3_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity3_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity3_mvaluesbottom
#define seglength mccelliptical_guide_gravity3_seglength
#define guideInfo mccelliptical_guide_gravity3_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity3_latestParticleCollision
#define Gx mccelliptical_guide_gravity3_Gx
#define Gy mccelliptical_guide_gravity3_Gy
#define Gz mccelliptical_guide_gravity3_Gz
#define Gx0 mccelliptical_guide_gravity3_Gx0
#define Gy0 mccelliptical_guide_gravity3_Gy0
#define Gz0 mccelliptical_guide_gravity3_Gz0
#define Circ mccelliptical_guide_gravity3_Circ
#define dynamicalSegLength mccelliptical_guide_gravity3_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity3=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity3_l;
MCNUM xwidth = mccelliptical_guide_gravity3_xwidth;
MCNUM yheight = mccelliptical_guide_gravity3_yheight;
MCNUM linxw = mccelliptical_guide_gravity3_linxw;
MCNUM loutxw = mccelliptical_guide_gravity3_loutxw;
MCNUM linyh = mccelliptical_guide_gravity3_linyh;
MCNUM loutyh = mccelliptical_guide_gravity3_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity3_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity3_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity3_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity3_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity3_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity3_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity3_dimensionsAt;
char* option = mccelliptical_guide_gravity3_option;
MCNUM R0 = mccelliptical_guide_gravity3_R0;
MCNUM Qc = mccelliptical_guide_gravity3_Qc;
MCNUM alpha = mccelliptical_guide_gravity3_alpha;
MCNUM m = mccelliptical_guide_gravity3_m;
MCNUM W = mccelliptical_guide_gravity3_W;
MCNUM alpharight = mccelliptical_guide_gravity3_alpharight;
MCNUM mright = mccelliptical_guide_gravity3_mright;
MCNUM alphaleft = mccelliptical_guide_gravity3_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity3_mleft;
MCNUM alphatop = mccelliptical_guide_gravity3_alphatop;
MCNUM mtop = mccelliptical_guide_gravity3_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity3_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity3_mbottom;
char* verbose = mccelliptical_guide_gravity3_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity3_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity3_curvature;
#line 1496 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{


	// Calculate the points need to draw approximation of the ellipses
	// defining the guide 

	// the number of lines used to draw one side of the guide
	int ApproximationMirrors = 500; 

	// The start of the guide
	double zvalue=0;

	// The the different in z between point used to draw the ellipse
	double zdelta = guideInfo.Length/(1.0*ApproximationMirrors);

	// The vector used to store the points defining the lines
	double xplus[ApproximationMirrors+1];
	double xminus[ApproximationMirrors+1];
	double yplus[ApproximationMirrors+1];
	double yminus[ApproximationMirrors+1];

	// Temperary values for the loop
	double tempx;
	double tempy;

	/*
		Calculate the second coordinates to the points on the ellipse with z_i 
		as the first coordinate. We transform the point to the coordinate system 
		there the ellipse is a unit circle. And use the defination of this circle 
		to find second coordinate (x^2+z^2 = 1)
	*/

	/////////////////////////////////////////////////////////
	double Length;
	double entranceHorizontalWidth;
	double entranceVerticalWidth;	
	
	// ellipses infomation
	double ellipseMajorAxis[4],		ellipseMinorAxis[4];
	double ellipseMajorOffset[4],	ellipseMinorOffset[4];

	enum Side {RightSide,TopSide,LeftSide,BottomSide,None};
	/////////////////////////////////////////////////////////

	int i = 0;
	double tempz = 0;
	for(i=0;i<ApproximationMirrors+1;i++){
		
		tempx = sqrt(
			guideInfo.ellipseMinorAxis[RightSide]
			*guideInfo.ellipseMinorAxis[RightSide]
			-(	guideInfo.ellipseMinorAxis[RightSide]
				*guideInfo.ellipseMinorAxis[RightSide] )
			/(	guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
		);

		xplus[i] =  tempx + guideInfo.ellipseMinorOffset[RightSide];
		xminus[i]= -tempx + guideInfo.ellipseMinorOffset[RightSide];
		
		tempy = sqrt(
			guideInfo.ellipseMinorAxis[TopSide]
			*guideInfo.ellipseMinorAxis[TopSide]
			-(	guideInfo.ellipseMinorAxis[TopSide]
				*guideInfo.ellipseMinorAxis[TopSide] )
			/(	guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
		);
		
		yplus[i] =  tempy + guideInfo.ellipseMinorOffset[TopSide];
		yminus[i]= -tempy + guideInfo.ellipseMinorOffset[TopSide];
	}

	///// Draw lines

	// Drawing lines orthogonal with the z direction.
	// at both ends of the guide and at the boardest place at the guide

	// These may not give correct result if one of the ends are closed

	int j=0;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 	  , yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta,xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j],0 		, zvalue+j*zdelta);
			
	j=ApproximationMirrors;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j],zvalue+j*zdelta,xminus[j] , yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	    , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		 , zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j], 0 		 , zvalue+j*zdelta);

			  
			  
	// find boardest place on the guide and draw a band around the guide
	int m0;
	double boardestPlace = 0;
	int boardestPlaceNumber = 0;
	for(m0=0; m0<ApproximationMirrors; m0++){
		if( boardestPlace <= fabs(yplus[m0]) ){
			boardestPlace = fabs(yplus[m0]);
			boardestPlaceNumber = m0;
		}
	}
	j = boardestPlaceNumber;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta, xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta, xplus[j], 0 		 , zvalue+j*zdelta);


			  
	// Drawing lines parallel with the z direction
			  
	int k=0;
	for(k=0; k < ApproximationMirrors; k++){
		
		line( xplus[k], yplus[k], zvalue+k*zdelta,xplus[k+1], yplus[k+1], zvalue+(k+1)*zdelta);
		line( xminus[k],yplus[k], zvalue+k*zdelta,xminus[k+1],yplus[k+1], zvalue+(k+1)*zdelta);
		
		line( xplus[k], yminus[k],zvalue+k*zdelta,xplus[k+1], yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],yminus[k],zvalue+k*zdelta,xminus[k+1],yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],0 		, zvalue+k*zdelta, xminus[k+1],0 	, zvalue+(k+1)*zdelta);
		line( xplus[k], 0 		, zvalue+k*zdelta, xplus[k+1], 0 	, zvalue+(k+1)*zdelta);
		
		line( 0 	, yminus[k],  zvalue+k*zdelta, 0 	  	,yminus[k+1],zvalue+(k+1)*zdelta);
		line( 0 	,yplus[k]  , zvalue+k*zdelta , 0 		,yplus[k] , zvalue+(k+1)*zdelta);
		
		}

		if(guideInfo.EnclosingBoxOn){
			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );

			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );
			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
		}
}
#line 19884 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity3=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EndOfelement_3'. */
  SIG_MESSAGE("EndOfelement_3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EndOfelement_3");
#define mccompcurname  EndOfelement_3
#define mccompcurtype  Arm
#define mccompcurindex 8
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 19919 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Curved_guide_2'. */
  SIG_MESSAGE("Curved_guide_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Curved_guide_2");
#define mccompcurname  Curved_guide_2
#define mccompcurtype  Bender
#define mccompcurindex 9
#define bk mccCurved_guide_2_bk
#define mWin mccCurved_guide_2_mWin
{   /* Declarations of Curved_guide_2=Bender() SETTING parameters. */
MCNUM w = mccCurved_guide_2_w;
MCNUM h = mccCurved_guide_2_h;
MCNUM r = mccCurved_guide_2_r;
MCNUM Win = mccCurved_guide_2_Win;
MCNUM k = mccCurved_guide_2_k;
MCNUM d = mccCurved_guide_2_d;
MCNUM l = mccCurved_guide_2_l;
MCNUM R0a = mccCurved_guide_2_R0a;
MCNUM Qca = mccCurved_guide_2_Qca;
MCNUM alphaa = mccCurved_guide_2_alphaa;
MCNUM ma = mccCurved_guide_2_ma;
MCNUM Wa = mccCurved_guide_2_Wa;
MCNUM R0i = mccCurved_guide_2_R0i;
MCNUM Qci = mccCurved_guide_2_Qci;
MCNUM alphai = mccCurved_guide_2_alphai;
MCNUM mi = mccCurved_guide_2_mi;
MCNUM Wi = mccCurved_guide_2_Wi;
MCNUM R0s = mccCurved_guide_2_R0s;
MCNUM Qcs = mccCurved_guide_2_Qcs;
MCNUM alphas = mccCurved_guide_2_alphas;
MCNUM ms = mccCurved_guide_2_ms;
MCNUM Ws = mccCurved_guide_2_Ws;
#line 269 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Bender.comp"
{
  int i;
  double w1c, w2c, h1, h2, L, w1, w2;

  w1c = (w + d)/(double)k;
  w2c = w1c; h1 = h; h2 = h;
  L = r*mWin; w1 = w; w2 = w;

  
  for(i = 0; i < k; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)L,
              i*w2c - w2/2.0,  h2/2.0, (double)L,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)L,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)L,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)L, w2/2.0, -h2/2.0, (double)L);
}
#line 19983 "./FullInstrument_v4.7_frontend.c"
}   /* End of Curved_guide_2=Bender() SETTING parameter declarations. */
#undef mWin
#undef bk
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EndOfelement_2'. */
  SIG_MESSAGE("EndOfelement_2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EndOfelement_2");
#define mccompcurname  EndOfelement_2
#define mccompcurtype  Arm
#define mccompcurindex 10
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20005 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'elliptical_guide_gravity1'. */
  SIG_MESSAGE("elliptical_guide_gravity1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "elliptical_guide_gravity1");
#define mccompcurname  elliptical_guide_gravity1
#define mccompcurtype  Elliptic_guide_gravity
#define mccompcurindex 11
#define mvaluesright mccelliptical_guide_gravity1_mvaluesright
#define mvaluesleft mccelliptical_guide_gravity1_mvaluesleft
#define mvaluestop mccelliptical_guide_gravity1_mvaluestop
#define mvaluesbottom mccelliptical_guide_gravity1_mvaluesbottom
#define seglength mccelliptical_guide_gravity1_seglength
#define guideInfo mccelliptical_guide_gravity1_guideInfo
#define latestParticleCollision mccelliptical_guide_gravity1_latestParticleCollision
#define Gx mccelliptical_guide_gravity1_Gx
#define Gy mccelliptical_guide_gravity1_Gy
#define Gz mccelliptical_guide_gravity1_Gz
#define Gx0 mccelliptical_guide_gravity1_Gx0
#define Gy0 mccelliptical_guide_gravity1_Gy0
#define Gz0 mccelliptical_guide_gravity1_Gz0
#define Circ mccelliptical_guide_gravity1_Circ
#define dynamicalSegLength mccelliptical_guide_gravity1_dynamicalSegLength
{   /* Declarations of elliptical_guide_gravity1=Elliptic_guide_gravity() SETTING parameters. */
MCNUM l = mccelliptical_guide_gravity1_l;
MCNUM xwidth = mccelliptical_guide_gravity1_xwidth;
MCNUM yheight = mccelliptical_guide_gravity1_yheight;
MCNUM linxw = mccelliptical_guide_gravity1_linxw;
MCNUM loutxw = mccelliptical_guide_gravity1_loutxw;
MCNUM linyh = mccelliptical_guide_gravity1_linyh;
MCNUM loutyh = mccelliptical_guide_gravity1_loutyh;
MCNUM majorAxisxw = mccelliptical_guide_gravity1_majorAxisxw;
MCNUM minorAxisxw = mccelliptical_guide_gravity1_minorAxisxw;
MCNUM majorAxisyh = mccelliptical_guide_gravity1_majorAxisyh;
MCNUM minorAxisyh = mccelliptical_guide_gravity1_minorAxisyh;
MCNUM majorAxisoffsetxw = mccelliptical_guide_gravity1_majorAxisoffsetxw;
MCNUM majorAxisoffsetyh = mccelliptical_guide_gravity1_majorAxisoffsetyh;
char* dimensionsAt = mccelliptical_guide_gravity1_dimensionsAt;
char* option = mccelliptical_guide_gravity1_option;
MCNUM R0 = mccelliptical_guide_gravity1_R0;
MCNUM Qc = mccelliptical_guide_gravity1_Qc;
MCNUM alpha = mccelliptical_guide_gravity1_alpha;
MCNUM m = mccelliptical_guide_gravity1_m;
MCNUM W = mccelliptical_guide_gravity1_W;
MCNUM alpharight = mccelliptical_guide_gravity1_alpharight;
MCNUM mright = mccelliptical_guide_gravity1_mright;
MCNUM alphaleft = mccelliptical_guide_gravity1_alphaleft;
MCNUM mleft = mccelliptical_guide_gravity1_mleft;
MCNUM alphatop = mccelliptical_guide_gravity1_alphatop;
MCNUM mtop = mccelliptical_guide_gravity1_mtop;
MCNUM alphabottom = mccelliptical_guide_gravity1_alphabottom;
MCNUM mbottom = mccelliptical_guide_gravity1_mbottom;
char* verbose = mccelliptical_guide_gravity1_verbose;
MCNUM enableGravity = mccelliptical_guide_gravity1_enableGravity;
MCNUM curvature = mccelliptical_guide_gravity1_curvature;
#line 1496 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Elliptic_guide_gravity.comp"
{


	// Calculate the points need to draw approximation of the ellipses
	// defining the guide 

	// the number of lines used to draw one side of the guide
	int ApproximationMirrors = 500; 

	// The start of the guide
	double zvalue=0;

	// The the different in z between point used to draw the ellipse
	double zdelta = guideInfo.Length/(1.0*ApproximationMirrors);

	// The vector used to store the points defining the lines
	double xplus[ApproximationMirrors+1];
	double xminus[ApproximationMirrors+1];
	double yplus[ApproximationMirrors+1];
	double yminus[ApproximationMirrors+1];

	// Temperary values for the loop
	double tempx;
	double tempy;

	/*
		Calculate the second coordinates to the points on the ellipse with z_i 
		as the first coordinate. We transform the point to the coordinate system 
		there the ellipse is a unit circle. And use the defination of this circle 
		to find second coordinate (x^2+z^2 = 1)
	*/

	/////////////////////////////////////////////////////////
	double Length;
	double entranceHorizontalWidth;
	double entranceVerticalWidth;	
	
	// ellipses infomation
	double ellipseMajorAxis[4],		ellipseMinorAxis[4];
	double ellipseMajorOffset[4],	ellipseMinorOffset[4];

	enum Side {RightSide,TopSide,LeftSide,BottomSide,None};
	/////////////////////////////////////////////////////////

	int i = 0;
	double tempz = 0;
	for(i=0;i<ApproximationMirrors+1;i++){
		
		tempx = sqrt(
			guideInfo.ellipseMinorAxis[RightSide]
			*guideInfo.ellipseMinorAxis[RightSide]
			-(	guideInfo.ellipseMinorAxis[RightSide]
				*guideInfo.ellipseMinorAxis[RightSide] )
			/(	guideInfo.ellipseMajorAxis[RightSide]
				*guideInfo.ellipseMajorAxis[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[RightSide] )
		);

		xplus[i] =  tempx + guideInfo.ellipseMinorOffset[RightSide];
		xminus[i]= -tempx + guideInfo.ellipseMinorOffset[RightSide];
		
		tempy = sqrt(
			guideInfo.ellipseMinorAxis[TopSide]
			*guideInfo.ellipseMinorAxis[TopSide]
			-(	guideInfo.ellipseMinorAxis[TopSide]
				*guideInfo.ellipseMinorAxis[TopSide] )
			/(	guideInfo.ellipseMajorAxis[TopSide]
				*guideInfo.ellipseMajorAxis[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
			*( zvalue+zdelta*i-guideInfo.ellipseMajorOffset[TopSide] )
		);
		
		yplus[i] =  tempy + guideInfo.ellipseMinorOffset[TopSide];
		yminus[i]= -tempy + guideInfo.ellipseMinorOffset[TopSide];
	}

	///// Draw lines

	// Drawing lines orthogonal with the z direction.
	// at both ends of the guide and at the boardest place at the guide

	// These may not give correct result if one of the ends are closed

	int j=0;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 	  , yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta,xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j],0 		, zvalue+j*zdelta);
			
	j=ApproximationMirrors;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j],zvalue+j*zdelta,xminus[j] , yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	    , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		 , zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta,xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta,xplus[j], 0 		 , zvalue+j*zdelta);

			  
			  
	// find boardest place on the guide and draw a band around the guide
	int m0;
	double boardestPlace = 0;
	int boardestPlaceNumber = 0;
	for(m0=0; m0<ApproximationMirrors; m0++){
		if( boardestPlace <= fabs(yplus[m0]) ){
			boardestPlace = fabs(yplus[m0]);
			boardestPlaceNumber = m0;
		}
	}
	j = boardestPlaceNumber;

	line( xplus[j], yplus[j], zvalue+j*zdelta,0 		, yplus[j], zvalue+j*zdelta);
	line( 0 	   , yplus[j], zvalue+j*zdelta,xminus[j], yplus[j], zvalue+j*zdelta);
	line( xminus[j], yminus[j], zvalue+j*zdelta,0 	   , yminus[j], zvalue+j*zdelta);
	line( 0, 		yminus[j], zvalue+j*zdelta, xplus[j], yminus[j], zvalue+j*zdelta);
	line( xminus[j],yplus[j], zvalue+j*zdelta, xminus[j],0 		, zvalue+j*zdelta);
	line( xminus[j],0 		 , zvalue+j*zdelta, xminus[j],yminus[j], zvalue+j*zdelta);
	line( xplus[j], 0, 		  zvalue+j*zdelta, xplus[j], yplus[j], zvalue+j*zdelta);
	line( xplus[j], yminus[j], zvalue+j*zdelta, xplus[j], 0 		 , zvalue+j*zdelta);


			  
	// Drawing lines parallel with the z direction
			  
	int k=0;
	for(k=0; k < ApproximationMirrors; k++){
		
		line( xplus[k], yplus[k], zvalue+k*zdelta,xplus[k+1], yplus[k+1], zvalue+(k+1)*zdelta);
		line( xminus[k],yplus[k], zvalue+k*zdelta,xminus[k+1],yplus[k+1], zvalue+(k+1)*zdelta);
		
		line( xplus[k], yminus[k],zvalue+k*zdelta,xplus[k+1], yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],yminus[k],zvalue+k*zdelta,xminus[k+1],yminus[k+1],zvalue+(k+1)*zdelta);
		
		line( xminus[k],0 		, zvalue+k*zdelta, xminus[k+1],0 	, zvalue+(k+1)*zdelta);
		line( xplus[k], 0 		, zvalue+k*zdelta, xplus[k+1], 0 	, zvalue+(k+1)*zdelta);
		
		line( 0 	, yminus[k],  zvalue+k*zdelta, 0 	  	,yminus[k+1],zvalue+(k+1)*zdelta);
		line( 0 	,yplus[k]  , zvalue+k*zdelta , 0 		,yplus[k] , zvalue+(k+1)*zdelta);
		
		}

		if(guideInfo.EnclosingBoxOn){
			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );

			dashed_line( guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],
							guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],10 );
			dashed_line( guideInfo.xArray[4],guideInfo.yArray[4],guideInfo.zArray[4],
							guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],10 );
			dashed_line( guideInfo.xArray[7],guideInfo.yArray[7],guideInfo.zArray[7],
							guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],10 );
			dashed_line( guideInfo.xArray[3],guideInfo.yArray[3],guideInfo.zArray[3],
							guideInfo.xArray[0],guideInfo.yArray[0],guideInfo.zArray[0],10 );

			dashed_line( guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],
							guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],10 );
			dashed_line( guideInfo.xArray[5],guideInfo.yArray[5],guideInfo.zArray[5],
							guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],10 );
			dashed_line( guideInfo.xArray[6],guideInfo.yArray[6],guideInfo.zArray[6],
							guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],10 );
			dashed_line( guideInfo.xArray[2],guideInfo.yArray[2],guideInfo.zArray[2],
							guideInfo.xArray[1],guideInfo.yArray[1],guideInfo.zArray[1],10 );
		}
}
#line 20253 "./FullInstrument_v4.7_frontend.c"
}   /* End of elliptical_guide_gravity1=Elliptic_guide_gravity() SETTING parameter declarations. */
#undef dynamicalSegLength
#undef Circ
#undef Gz0
#undef Gy0
#undef Gx0
#undef Gz
#undef Gy
#undef Gx
#undef latestParticleCollision
#undef guideInfo
#undef seglength
#undef mvaluesbottom
#undef mvaluestop
#undef mvaluesleft
#undef mvaluesright
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'EndOfelement_1'. */
  SIG_MESSAGE("EndOfelement_1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "EndOfelement_1");
#define mccompcurname  EndOfelement_1
#define mccompcurtype  Arm
#define mccompcurindex 12
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20288 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SlitGuideEnd'. */
  SIG_MESSAGE("SlitGuideEnd (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SlitGuideEnd");
#define mccompcurname  SlitGuideEnd
#define mccompcurtype  Slit
#define mccompcurindex 13
{   /* Declarations of SlitGuideEnd=Slit() SETTING parameters. */
MCNUM xmin = mccSlitGuideEnd_xmin;
MCNUM xmax = mccSlitGuideEnd_xmax;
MCNUM ymin = mccSlitGuideEnd_ymin;
MCNUM ymax = mccSlitGuideEnd_ymax;
MCNUM radius = mccSlitGuideEnd_radius;
MCNUM xwidth = mccSlitGuideEnd_xwidth;
MCNUM yheight = mccSlitGuideEnd_yheight;
#line 83 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
  
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
}
#line 20330 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitGuideEnd=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SlitMonitor'. */
  SIG_MESSAGE("SlitMonitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SlitMonitor");
#define mccompcurname  SlitMonitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 14
#define PSD_N mccSlitMonitor_PSD_N
#define PSD_p mccSlitMonitor_PSD_p
#define PSD_p2 mccSlitMonitor_PSD_p2
{   /* Declarations of SlitMonitor=PSD_monitor() SETTING parameters. */
int nx = mccSlitMonitor_nx;
int ny = mccSlitMonitor_ny;
char* filename = mccSlitMonitor_filename;
MCNUM xmin = mccSlitMonitor_xmin;
MCNUM xmax = mccSlitMonitor_xmax;
MCNUM ymin = mccSlitMonitor_ymin;
MCNUM ymax = mccSlitMonitor_ymax;
MCNUM xwidth = mccSlitMonitor_xwidth;
MCNUM yheight = mccSlitMonitor_yheight;
MCNUM restore_neutron = mccSlitMonitor_restore_neutron;
int nowritefile = mccSlitMonitor_nowritefile;
#line 132 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 20366 "./FullInstrument_v4.7_frontend.c"
}   /* End of SlitMonitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'MonoFLat'. */
  SIG_MESSAGE("MonoFLat (McDisplay)");
  printf("MCDISPLAY: component %s\n", "MonoFLat");
#define mccompcurname  MonoFLat
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 15
#define mos_rms_y mccMonoFLat_mos_rms_y
#define mos_rms_z mccMonoFLat_mos_rms_z
#define mos_rms_max mccMonoFLat_mos_rms_max
#define mono_Q mccMonoFLat_mono_Q
#define SlabWidth mccMonoFLat_SlabWidth
#define SlabHeight mccMonoFLat_SlabHeight
#define rTable mccMonoFLat_rTable
#define tTable mccMonoFLat_tTable
#define row mccMonoFLat_row
#define col mccMonoFLat_col
#define tiltH mccMonoFLat_tiltH
#define tiltV mccMonoFLat_tiltV
{   /* Declarations of MonoFLat=Monochromator_curved() SETTING parameters. */
char* reflect = mccMonoFLat_reflect;
char* transmit = mccMonoFLat_transmit;
MCNUM zwidth = mccMonoFLat_zwidth;
MCNUM yheight = mccMonoFLat_yheight;
MCNUM gap = mccMonoFLat_gap;
MCNUM NH = mccMonoFLat_NH;
MCNUM NV = mccMonoFLat_NV;
MCNUM mosaich = mccMonoFLat_mosaich;
MCNUM mosaicv = mccMonoFLat_mosaicv;
MCNUM r0 = mccMonoFLat_r0;
MCNUM t0 = mccMonoFLat_t0;
MCNUM Q = mccMonoFLat_Q;
MCNUM RV = mccMonoFLat_RV;
MCNUM RH = mccMonoFLat_RH;
MCNUM DM = mccMonoFLat_DM;
MCNUM mosaic = mccMonoFLat_mosaic;
MCNUM width = mccMonoFLat_width;
MCNUM height = mccMonoFLat_height;
MCNUM verbose = mccMonoFLat_verbose;
MCNUM order = mccMonoFLat_order;
#line 468 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  int ih;

  for(ih = 0; ih < NH; ih++)
  {
    int iv;
    for(iv = 0; iv < NV; iv++)
    {
      double zmin,zmax,ymin,ymax;
      double xt, yt;

      zmin = (SlabWidth+gap)*(ih-NH/2.0)+gap/2;
      zmax = zmin+SlabWidth;
      ymin = (SlabHeight+gap)*(iv-NV/2.0)+gap/2;
      ymax = ymin+SlabHeight;

      if (RH) xt = -(zmax*zmax - zmin*zmin)/RH/2;
      else    xt = 0;

      if (RV) yt = -(ymax*ymax - ymin*ymin)/RV/2;
      else    yt = 0;
      multiline(5, xt+yt, (double)ymin, (double)zmin,
                   xt-yt, (double)ymax, (double)zmin,
                  -xt-yt, (double)ymax, (double)zmax,
                  -xt+yt, (double)ymin, (double)zmax,
                   xt+yt, (double)ymin, (double)zmin);
     }
   }
}
#line 20444 "./FullInstrument_v4.7_frontend.c"
}   /* End of MonoFLat=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Mono'. */
  SIG_MESSAGE("Mono (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Mono");
#define mccompcurname  Mono
#define mccompcurtype  Monochromator_curved
#define mccompcurindex 16
#define mos_rms_y mccMono_mos_rms_y
#define mos_rms_z mccMono_mos_rms_z
#define mos_rms_max mccMono_mos_rms_max
#define mono_Q mccMono_mono_Q
#define SlabWidth mccMono_SlabWidth
#define SlabHeight mccMono_SlabHeight
#define rTable mccMono_rTable
#define tTable mccMono_tTable
#define row mccMono_row
#define col mccMono_col
#define tiltH mccMono_tiltH
#define tiltV mccMono_tiltV
{   /* Declarations of Mono=Monochromator_curved() SETTING parameters. */
char* reflect = mccMono_reflect;
char* transmit = mccMono_transmit;
MCNUM zwidth = mccMono_zwidth;
MCNUM yheight = mccMono_yheight;
MCNUM gap = mccMono_gap;
MCNUM NH = mccMono_NH;
MCNUM NV = mccMono_NV;
MCNUM mosaich = mccMono_mosaich;
MCNUM mosaicv = mccMono_mosaicv;
MCNUM r0 = mccMono_r0;
MCNUM t0 = mccMono_t0;
MCNUM Q = mccMono_Q;
MCNUM RV = mccMono_RV;
MCNUM RH = mccMono_RH;
MCNUM DM = mccMono_DM;
MCNUM mosaic = mccMono_mosaic;
MCNUM width = mccMono_width;
MCNUM height = mccMono_height;
MCNUM verbose = mccMono_verbose;
MCNUM order = mccMono_order;
#line 468 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Monochromator_curved.comp"
{
  int ih;

  for(ih = 0; ih < NH; ih++)
  {
    int iv;
    for(iv = 0; iv < NV; iv++)
    {
      double zmin,zmax,ymin,ymax;
      double xt, yt;

      zmin = (SlabWidth+gap)*(ih-NH/2.0)+gap/2;
      zmax = zmin+SlabWidth;
      ymin = (SlabHeight+gap)*(iv-NV/2.0)+gap/2;
      ymax = ymin+SlabHeight;

      if (RH) xt = -(zmax*zmax - zmin*zmin)/RH/2;
      else    xt = 0;

      if (RV) yt = -(ymax*ymax - ymin*ymin)/RV/2;
      else    yt = 0;
      multiline(5, xt+yt, (double)ymin, (double)zmin,
                   xt-yt, (double)ymax, (double)zmin,
                  -xt-yt, (double)ymax, (double)zmax,
                  -xt+yt, (double)ymin, (double)zmax,
                   xt+yt, (double)ymin, (double)zmin);
     }
   }
}
#line 20531 "./FullInstrument_v4.7_frontend.c"
}   /* End of Mono=Monochromator_curved() SETTING parameter declarations. */
#undef tiltV
#undef tiltH
#undef col
#undef row
#undef tTable
#undef rTable
#undef SlabHeight
#undef SlabWidth
#undef mono_Q
#undef mos_rms_max
#undef mos_rms_z
#undef mos_rms_y
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'A2'. */
  SIG_MESSAGE("A2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "A2");
#define mccompcurname  A2
#define mccompcurtype  Arm
#define mccompcurindex 17
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20563 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit'. */
  SIG_MESSAGE("slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit");
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 18
{   /* Declarations of slit=Slit() SETTING parameters. */
MCNUM xmin = mccslit_xmin;
MCNUM xmax = mccslit_xmax;
MCNUM ymin = mccslit_ymin;
MCNUM ymax = mccslit_ymax;
MCNUM radius = mccslit_radius;
MCNUM xwidth = mccslit_xwidth;
MCNUM yheight = mccslit_yheight;
#line 83 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
  
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
}
#line 20605 "./FullInstrument_v4.7_frontend.c"
}   /* End of slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'SamplePos'. */
  SIG_MESSAGE("SamplePos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "SamplePos");
#define mccompcurname  SamplePos
#define mccompcurtype  Arm
#define mccompcurindex 19
#line 40 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 20625 "./FullInstrument_v4.7_frontend.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_sample_pos'. */
  SIG_MESSAGE("psd_monitor_sample_pos (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_sample_pos");
#define mccompcurname  psd_monitor_sample_pos
#define mccompcurtype  PSD_monitor
#define mccompcurindex 20
#define PSD_N mccpsd_monitor_sample_pos_PSD_N
#define PSD_p mccpsd_monitor_sample_pos_PSD_p
#define PSD_p2 mccpsd_monitor_sample_pos_PSD_p2
{   /* Declarations of psd_monitor_sample_pos=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_sample_pos_nx;
int ny = mccpsd_monitor_sample_pos_ny;
char* filename = mccpsd_monitor_sample_pos_filename;
MCNUM xmin = mccpsd_monitor_sample_pos_xmin;
MCNUM xmax = mccpsd_monitor_sample_pos_xmax;
MCNUM ymin = mccpsd_monitor_sample_pos_ymin;
MCNUM ymax = mccpsd_monitor_sample_pos_ymax;
MCNUM xwidth = mccpsd_monitor_sample_pos_xwidth;
MCNUM yheight = mccpsd_monitor_sample_pos_yheight;
MCNUM restore_neutron = mccpsd_monitor_sample_pos_restore_neutron;
int nowritefile = mccpsd_monitor_sample_pos_nowritefile;
#line 132 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 20660 "./FullInstrument_v4.7_frontend.c"
}   /* End of psd_monitor_sample_pos=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'e_monitor_sample_ps'. */
  SIG_MESSAGE("e_monitor_sample_ps (McDisplay)");
  printf("MCDISPLAY: component %s\n", "e_monitor_sample_ps");
#define mccompcurname  e_monitor_sample_ps
#define mccompcurtype  E_monitor
#define mccompcurindex 21
#define nE mcce_monitor_sample_ps_nE
#define E_N mcce_monitor_sample_ps_E_N
#define E_p mcce_monitor_sample_ps_E_p
#define E_p2 mcce_monitor_sample_ps_E_p2
#define S_p mcce_monitor_sample_ps_S_p
#define S_pE mcce_monitor_sample_ps_S_pE
#define S_pE2 mcce_monitor_sample_ps_S_pE2
{   /* Declarations of e_monitor_sample_ps=E_monitor() SETTING parameters. */
char* filename = mcce_monitor_sample_ps_filename;
MCNUM xmin = mcce_monitor_sample_ps_xmin;
MCNUM xmax = mcce_monitor_sample_ps_xmax;
MCNUM ymin = mcce_monitor_sample_ps_ymin;
MCNUM ymax = mcce_monitor_sample_ps_ymax;
MCNUM xwidth = mcce_monitor_sample_ps_xwidth;
MCNUM yheight = mcce_monitor_sample_ps_yheight;
MCNUM Emin = mcce_monitor_sample_ps_Emin;
MCNUM Emax = mcce_monitor_sample_ps_Emax;
MCNUM restore_neutron = mcce_monitor_sample_ps_restore_neutron;
int nowritefile = mcce_monitor_sample_ps_nowritefile;
#line 132 "/usr/share/mcstas/2.7.1/tools/Python/mcrun/../mccodelib/../../../monitors/E_monitor.comp"
{
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
}
#line 20703 "./FullInstrument_v4.7_frontend.c"
}   /* End of e_monitor_sample_ps=E_monitor() SETTING parameter declarations. */
#undef S_pE2
#undef S_pE
#undef S_p
#undef E_p2
#undef E_p
#undef E_N
#undef nE
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
#undef cylinder
#undef sphere
/* end of generated C code ./FullInstrument_v4.7_frontend.c */
