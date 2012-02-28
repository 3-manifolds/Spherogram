/* morwen_repack.c 					8/95 */

/* This is a specialized repack engine to be called by Morwen 
Thistlethwaite's knot package. Intended only for small hyperbolic 
maximal packings. Looks for data in "/tmp/knot-comb.p" and 
returns a euclidean packing in "/tmp/knot-packing.p".

Compile with 

	cc -o mr morwen_repack.c -g -lm 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 5001		/* max number of vertices of packing */
#define ITERATIONS 1000000	/* max passes through the repack routine */
#define NUM_PETALS 60		/* max number of neighbors of any circle */
#define MAX_COMPONENTS 1	/* number of components of complex should be simply connected here */
#define TOLER .000001		/* threshold for repack computations */
#define OKERR .00000001		/* Roundoff for zero */
#define float double

/* Debug output control */
int OFF = 0;
int SUMMARY = 1;
int INTERMEDIATE = 2;
int DETAIL = 3;
extern int printf_bool;

/* geometry control */
extern int DRAW_IN_HYPERBOLIC_DISC;

/* ---- standard CirclePack data structures ------- */

typedef struct {float re,im;}  complex;

typedef struct
 {
	int vert[3];		/* ordered triple of vertices */
	int index_flag;		/* which to draw first, 0,1, or 2 */
	int next_face;		/* next face in drawing order. */
	int rwb_flag;		/* red/white/blue faces:0=white,1=red,2=blue see info file for description */
	int next_red;		/* next red face */
	int plot_flag;		/* useful flag */
	int color;		/* color index */
	int mark;		/* mark */
  } f_data;

struct Vertlist
 {
	int v;			/* vertex number */
	struct Vertlist *next;
 };

struct K_data
 {
	int num;                /* number of faces this node is in.*/
				/* 'flower' has entries from 0 to num */
	int bdry_flag; 		/* true for boundary nodes */
	int plot_flag;		/* true if node has been plotted - used when circles are being plotted */
	int flower[NUM_PETALS+1];	/* list of nodes in flower of this node, in positive orientation, last = first if flower is complete */
	int color;		/* gives color value, 0=background */
	int mark;		/* mark */
 } packK[MAX_SIZE];

struct R_data
 {
	complex center;		/* center as complex number */
	float rad;		/* radius of circle for this node. hyperbolic case: called 's_radius' for special. exp(-r) is stored for
				       finite radii; for infinite hyp radius, we store the negative of the eucl radius (can then use eucl data for plotting). */
	float curv;		/* angle sum at this vertex. */
	float aim;		/* desired curvature at this vertex */
 } packR[MAX_SIZE];

/* ==================== structure for pack data ====================*/

struct p_data
 {
	int status;		/* 0 if pack empty, 1 otherwise */
	char fhead[64];		/* pack header <= 40 char*/
	int hes;		/* curvature of geometry,-1=hyp,0=eucl,1=sph */
	struct K_data *packK_ptr;	/* pointer to corres complex data */
	struct R_data *packR_ptr;	/* pointer to corres rad data */
	int nodecount;		/* number of (finite) nodes */
	int facecount;		/* number of faces = 1+edges-nodes */
	int intnode;		/* number of interior nodes */
	int num_bdry_comp;	/* number of bdry components */
	int num_int_comp;	/* number of components of interior verts */
	int bdry_starts[MAX_COMPONENTS+1];	
				/* indices for nodes marking bdry components */
	int int_starts[MAX_COMPONENTS+1];	
				/* indices of nodes marking int components */
	int first_face;		/* index of first face to plot */
	int first_red_face;	/* index of "red" chain of faces */
	int beta;		/* index of beta node (on boundary) */
	int alpha;		/* index of alpha node (origin) */
	int gamma;		/* index of node to be plotted on y>0 axis */
	int euler;		/* Euler characteristic */
	int genus;		/* genus of complex. Euler+#bdry=2-2g */
	int active_node;	/* currently active_node */
	int locks;		/* locks placed by remote processes;
					bitwise flags */
	f_data *faces;		/* pointer to face data */
 } packdata;

struct stface
 {
	int face;
	struct stface *prev;
	struct stface *next;
 };

/* ----- user specified variables set before compiling -------- */

FILE *fp;
char *outputfile;
int iterates=2;			/* parameter of current packing routines */
int *face_org;	/* temp list of faces by vertex */

/*extern*/ float cAbs(),h_comp_cos(),h_cos_overlap(),h_radcalc(),e_cos_overlap();
/*extern*/ complex cdiv(),cmult(),cconj(),csub(),cadd(),mob_norm_inv();
/*extern*/ complex mob_norm(),mob_trans();

/* -------- main routine ----------- */

circle_pack(infile, outfile)
char *infile, *outfile;

{
	int i,nextcount;

/* -------- read data phase ----------------- */

	packdata.packK_ptr=packK;
	packdata.packR_ptr=packR;
    /*    if ((fp=fopen("/tmp/knot-comb.p","r"))==NULL || !readpack()) */
        if ((fp=fopen(infile,"r"))==NULL || !readpack())
	 {
		fprintf(stderr,
		   "Failed to find or properly read '/tmp/knot-comb.p'.");
		exit(1);
	 }
	fclose(fp);
//	if (!complex_count(&packdata) || !facedraworder(&packdata,NULL))
	if (!complex_count(&packdata))
	 {
		fprintf(stderr,"Failed due to complex_count.");
		exit(2);
	 }

	if (!facedraworder(&packdata,NULL))
	 {
		fprintf(stderr,"Failed due to facedraworder.");
		exit(2);
	 }
/* ---------- repack phase ------------------ */

if (printf_bool == INTERMEDIATE)
{
	printf("initial vertex radii:\n");
    for (i=1;i<=packdata.nodecount;i++)
		printf("  vertex%d radius = %.6f\n",i,packR[i].rad);
}

	nextcount=h_riffle(&packdata,ITERATIONS);

if (printf_bool >= INTERMEDIATE)
{
	printf("final vertex radii:\n");
    for (i=1;i<=packdata.nodecount;i++)
		printf("  vertex%d radius = %.6f\n",i,packR[i].rad);
}
	if (nextcount>=ITERATIONS)
	 {
		fprintf(stderr,"Warning: repacking may be incomplete.");
	 }
	comp_pack_centers(&packdata);

	
/* ---------- convert data to euclidean ------ */

    if (!DRAW_IN_HYPERBOLIC_DISC)
    {
		packdata.hes=0;
		for (i=1;i<=packdata.nodecount;i++)
		{
			h_to_e_data(packR[i].center,packR[i].rad,
				&(packR[i].center),&(packR[i].rad)); 
		}
	}

/* ----------- write output phase ------------ */
        outputfile = outfile;
        
        writepack();

	return 0; /* success */
} /* end of main */

/* ================= input/output routines (minimal) =================== */

int
readpack() /* read pack "knot-comb.p". Return 1 if okay. */
{
	int i,k,count,j;
	char command[64];
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=packdata.packK_ptr;
	pR_ptr=packdata.packR_ptr;
	while ( (fscanf(fp,"%s",command)!=EOF)
		&& (strcmp(command,"END")!=0) )
	 {
		if (strcmp(command,"NODECOUNT:")==0) /* new packing */
		 {
			packdata.hes=-1; /* default */
			if (fscanf(fp,"%d",&count)!=1
				|| count<3 || count>MAX_SIZE)
				return -1;
			packdata.nodecount=count;
	 	 }
		if (strcmp(command,"ALPHA/BETA/GAMMA:")==0)
		 		fscanf(fp,"%d %d %d",&(packdata.alpha),
				&(packdata.beta),&(packdata.gamma));
		else if (strcmp(command,"FLOWERS:")==0)
		 {
		      for (i=1;i<=count;i++)
		      {
				  if (!fscanf(fp,"%d",&k)) // read the node index k
		              return 0;
		             
				  if ( !fscanf(fp,"%d",&pK_ptr[k].num) || pK_ptr[k].num>NUM_PETALS )  // read the number of petals around node k
				      return 0;
				      
				      
		          for (j=0;j<=pK_ptr[k].num;j++)
			      {
					   if (!fscanf(fp,"%d",&pK_ptr[k].flower[j])) // read the flower
						  return 0;
				  }
		      }
		 }
	} /* end of while */
	for (i=1;i<=count;i++)
	 {	
		pK_ptr[i].bdry_flag=0;
		pR_ptr[i].rad=0.1;
		pR_ptr[i].aim=6.28318530718; // 2\pi the target angle sum at the vertex.

//if (i == 1 || i == 5 || i == 10 || i==18)
//if (i==29 || i==30)
//	pR_ptr[i].aim *= 0.5;


		if (pK_ptr[i].flower[0]!=pK_ptr[i].flower[pK_ptr[i].num])
		 {
			pK_ptr[i].bdry_flag=1;
			pR_ptr[i].rad=-0.1;
			pR_ptr[i].aim=-0.1;
		 }
	 }
	for (i=1;i<=count;i++) 
		pR_ptr[i].center.re=pR_ptr[i].center.im=0.0;
	for (i=1;i<=count;i++)
		h_anglesum(&packdata,i,pR_ptr[i].rad, &(pR_ptr[i].curv));
	return (1);
} /* readpack */


int
writepack() /* write pack */
{
	int i,count,n;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=packdata.packK_ptr;
	pR_ptr=packdata.packR_ptr;
	count=packdata.nodecount;

	/*=====================insert by MBT=====================*/
	{
	FILE *fp1;
	fp1=fopen(outputfile , "w");
	
	int AB_type123_count = count;
    for (i=1;i<=count;i++)
		AB_type123_count -= pK_ptr[i].bdry_flag;
//	fprintf(fp1,"%d\n",count);
	fprintf(fp1,"%d\n",AB_type123_count);
	 
	
    for (i=1;i<=AB_type123_count;i++)
	{
		fprintf(fp1," % .6f % .6f  ",pR_ptr[i].center.re, pR_ptr[i].center.im);
		fprintf(fp1,"\n");
    }

	/* AB: write the radii as well, so metapost can draw circles */
    for (i=1;i<=AB_type123_count;i++)
		fprintf(fp1," %.6f\n",pR_ptr[i].rad);

    fclose(fp1);
    }
    /*=======================================================*/
	
	
	return 1;
} /* writepack */

/*=================== repacking routine =======================*/

int
h_riffle(p,passes) /* adjust s_radii to meet curvature targets in 'aim'.
aim<0 means that radius is fixed - no adjustments.
Bdry radii adjusted only if their aim >= 0. */
int passes;
struct p_data *p;
{
	int i,j,aimnum=0,*index,count=0,dummy;
	float recip,accum,verr,err,cut;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	index=(int *)malloc(5*(MAX_SIZE+1)*sizeof(int));
	pK_ptr=p->packK_ptr;pR_ptr=p->packR_ptr;
	accum=0;
	
	/* Since readpack sets the target angle sum aim to be 2\pi for interior vertices and -0.1
	   for boundary vertices, the following for-loop never matches the first condition, and always
	   matches the second for interior vertices.
	*/
	
	for (i=1;i<=p->nodecount;i++)
	{
		/* since boundary radii are initialized to -0.1, the first part of
		   the following clause is never executed
		*/
		if (pK_ptr[i].bdry_flag && pR_ptr[i].aim>=0  && pR_ptr[i].aim<.001)
		{

if (printf_bool)
printf("node %2d is a boundary node with small positive target angle sum, setting angle sum to -0.2\n",i);
			pR_ptr[i].rad=(-.2);
		}
		else if (pR_ptr[i].aim>0)
		{
			/* since interior radii are initialized to 0.1 the 
			   following clause is never executed
			*/
			if (pR_ptr[i].rad<=0 && pR_ptr[i].aim>.00001) 
			{
if (printf_bool)
printf("node %2d is an interior node with non-positive radius and target angle sum > 0.00001, setting radius to 0.01\n",i);
				pR_ptr[i].rad=.01;
			}
			
			/* index [k] = k+1 for all k, since the interior vertices appear first in the trinagulation data 
			   the pack data is indexed from 1 not zero because the vertex numbers in the triangulation data
			   are numbered from 1.
			*/
			index[aimnum]=i;
if (printf_bool >= DETAIL)
printf("index[%d]=%d\n",aimnum,i);			
			aimnum++;
			err=pR_ptr[i].curv-pR_ptr[i].aim;
			accum += (err<0) ? (-err) : err;
		}
	}
	
	if (aimnum==0) {free(index);return count;}
	recip=.333333/aimnum;
	
	/* we have cut = accnum/(3*num_interior_vertices) where accnum is the aggregate absolute delta
	   between the initial angle sum and the target angle sum.  The initial angle sum was calculated
	   by readpack from the initial radii.
	*/
	cut=accum*recip;

	while (cut >TOLER && count<passes)
	{

if (printf_bool)
printf("starting iteration with angle sum error threshold = %f\n",cut);

		for (j=0;j<aimnum;j++)
	    {
			i=index[j];
if (printf_bool)
printf("  node %2d ",i);
			h_anglesum(p,i,pR_ptr[i].rad,&pR_ptr[i].curv);
			verr=pR_ptr[i].curv-pR_ptr[i].aim;
if (printf_bool)
printf("current angle sum = %f, sum-2\\pi error = %f: ",pR_ptr[i].curv,verr);

			if (fabs(verr)>cut)
			{
			   pR_ptr[i].rad=h_radcalc(p,i,	pR_ptr[i].rad,pR_ptr[i].aim,&dummy);
if (printf_bool)
printf("radius adjusted to %f\n",pR_ptr[i].rad);
			   count++;
			}
			else
			{
if (printf_bool)
printf("radius not adjusted\n");
			}
		}
		
		accum=0;
		for (j=0;j<aimnum;j++)
		{
			i=index[j];
			err=pR_ptr[i].curv-pR_ptr[i].aim;
			accum += (err<0) ? (-err) : err;
		}
		
		cut=accum*recip;
	} /* end of while */
	
	free(index);
	return count;
} /* h_riffle */


/*==================== hyp and geometry routines ===========*/

float
h_radcalc(p,i,s,aim)  /* returns best radius at vertex i to
achieve anglesum aim. Use neighboring s_radii from pack and starting 
s_radius s (first guess). This routine is a mixture of binary search and
then Newton's method. `iterates' gives the max number of iterations to
use in trying for solution. */
struct p_data *p;
int i;
float s,aim;
{
	int n=0;
	float err,bestrad, bestcurv, curvderiv, tryrad, trycur;

	bestrad=s;
	
if (printf_bool >= DETAIL)
printf("h_radcalc called with current radius = %f, target = %f\n", bestrad, aim);

	if (aim<0.00001 && packK[i].bdry_flag) 
	   return (-.2);
	   
	if (s<=0 && aim>0.00001) 
	    bestrad=.5;
	    
	/* evaluate the  current angle sum in bestcurv and it's derivative in curvderiv */
	h_curvcalc(p,i,bestrad,&bestcurv,&curvderiv); 
	
if (printf_bool >= DETAIL)
printf("  current angle sum = %f, initial f(r) = %f\n", bestcurv, bestcurv-aim);

    /* The radius calculation is based on the function f(r) = \theta(r;{r_i}) - 2\pi, the delta between
       the angle sum and 2\pi
       
       the routine is looking for a radius that minimizes f (similar to step 3 of section 3.3 but using 
       \theta not \hat\theta)
       
       while |f(current_radius)| > TOLER update the current radius by doing Newton-Raphson then 
       some binary searching.
    */
	while ( (bestcurv>(aim+TOLER) || bestcurv<(aim-TOLER))
		 && (n < iterates))
	 {
		 
		/* do Newton-Raphson on f(r) = \theta(r;{r_i}) - 2\pi starting from r = current radius */
		tryrad=bestrad-(bestcurv-aim)/curvderiv;
		
if (printf_bool >= DETAIL)
printf("  Newton radius = %f\n", tryrad);

		if (tryrad<bestrad*0.5) 
		    tryrad=bestrad*0.5;
		    
		if (tryrad>(1+bestrad)*0.5) 
		    tryrad=(1+bestrad)*0.5;
		
if (printf_bool >= DETAIL)
printf("  initial trial radius  = adjusted Newton radius = %f\n", tryrad);

		/* evaluate angle sum with Newton Raphson radius in trycur */
		h_curvcalc(p,i,tryrad,&trycur,&curvderiv);
		
		/* set err to the absolute delta between the current angle sum and 2\pi 
		   that is err = |f(current_radius)|
		*/
		err=bestcurv-aim;
		
		if (err<0) 
		    err=(-err);
		
if (printf_bool >= DETAIL)
printf("  err threshold = f(current_radius) = f(%f) = %f\n", bestrad, bestcurv-aim);

		/* while |f(trycur)| > err, set the trial radius to the value halfway between
		   the trial value and the current radius
		   
		   DOESN'T THIS ALWAYS GO THE WRONG WAY?  In any case, the Newton radius is closer to the
		   root of f than the initial value, so surely we are never going to get |f(trycur)| > err,
		   so this code has no effect.
		*/
		while ( ((trycur-aim)>err || (trycur-aim)<(-err))
			 && (n<iterates))
		 {
			n++;
			tryrad=(tryrad+bestrad)*0.5;
			h_curvcalc(p,i,tryrad,&trycur,&curvderiv);
if (printf_bool >= DETAIL)
printf("    setting trial radius to %f gives f(trial_radius) = %f\n", tryrad, trycur-aim);
		 }
		n++;
		bestrad=tryrad;
		h_curvcalc(p,i,bestrad, &bestcurv,&curvderiv);

if (printf_bool >= DETAIL)
printf("  set current_radius to %f, with angle sum %f and f(current_radius) = %f\n", bestrad, bestcurv, bestcurv-aim);
	 }
	return bestrad;
} /* h_radcalc */

h_curvcalc(p,i,s,c,d)  /* comp angle sum c and its derivative d (with
respect to the s_radius, for interior nodes in hyp setting. Deriv
depends on whether some of radii are infinite */
struct p_data *p;
int i;
float s,*c,*d;
{
	int k;
	float s1,s2,a,b,cc,e,ss,cs;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=p->packK_ptr;pR_ptr=p->packR_ptr;
	*c=0;*d=0;
	if (s<=0) return; /* infinite radius at vertex of interest */
	s2=pR_ptr[pK_ptr[i].flower[0]].rad;
	for (k=1;k<=pK_ptr[i].num;k++)
	 {
		s1=s2;
		s2=pR_ptr[pK_ptr[i].flower[k]].rad;
		cs=h_comp_cos(s,s1,s2);
		*c+=acos(cs);
		if ((s1>0) && (s2>0))
		 {
			cc=(-s1*s1-s2*s2);
			b=s1*s1*s2*s2;
			a=(-2-cc-2*b);
			e=cc-a;
			ss=s*s;
			*d+=2*e*(1-b*ss*ss)/((1+cc*ss+b*ss*ss)*
				sqrt( e*(2+(cc+a)*ss+2*b*ss*ss) ));
		 }
		else if ((s1<=0) && (s2<=0)) *d+=2/sqrt(1-s*s);
		else if (s1<=0) *d+=sqrt(s2/s)/(1+s2*s);
		else *d+=sqrt(s1/s)/(1+s1*s);
	 }
} /*h_curvcalc */

int
h_anglesum(p,i,s,c) /* compute ang sum. */
int i;
float s,*c;
struct p_data *p;
{
	int k,m,j2,j1;
	float s1,s2,o1,o2,o3;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=p->packK_ptr;
	pR_ptr=p->packR_ptr;
	*c=0;
	if (s<=0) return; /* infinite radius at vertex of interest */
	j2=pK_ptr[i].flower[0];
	s2=pR_ptr[j2].rad;
	for (k=1;k<=pK_ptr[i].num;k++)
	 {
		s1=s2;
		s2=pR_ptr[pK_ptr[i].flower[k]].rad;
		*c+=acos(h_comp_cos(s,s1,s2));
	 }
} /* h_anglesum */

float
h_comp_cos(s1,s2,s3) /* given ordered triple of s_radii, compute
the cosine of the angle at first circle in triangle formed by
mutually tangent circles. */
float s1,s2,s3;
{
	float s1s,s2s,s3s,ans;

	if (s1<=0) return (1.0);
	s1s=s1*s1; 
	if ((s2<=0) && (s3<=0)) return (1-2*s1s);
	s3s=s3*s3;
	if (s2<=0) return ((1+s1s*(s3s-2))/(1-s1s*s3s));
	s2s=s2*s2;
 	if (s3<=0) return ((1+s1s*(s2s-2))/(1-s1s*s2s));
 	
	ans=((1+s1s*s2s)*(1+s1s*s3s)-2*s1s*(1+s2s*s3s))/
		((1-s1s*s2s)*(1-s1s*s3s));
	if (ans>1) return 1;
	if (ans<-1) return -1;
	return (ans);
} /* h_comp_cos */

int
h_to_e_data(h_center,s_rad,e_center,e_rad) /* converts circle data
from hyp to eucl. */
complex h_center;
float s_rad;
complex *e_center;
float *e_rad;
{
	float a,b,d,k,aec,ahc,g;

	if (s_rad<=0) /* infinite hyp radius */
	 {
		a=1+s_rad;
		e_center->re=h_center.re*a;
		e_center->im=h_center.im*a;
		*e_rad=(-s_rad); /* assumes -s_rad is meaningful */
		return 1;
	 }
	ahc=cAbs(h_center);
	g=((1+ahc)/(1-ahc));
	a=g/s_rad;
	d=(a-1)/(a+1);
	b=g*s_rad;
	k=(b-1)/(b+1);
	*e_rad = (d-k)*0.5;
	aec=(d+k)*0.5;
	if (ahc<=OKERR)
	 {
		e_center->re=0;
		e_center->im=0;
	 }
	else
	 {
		b=aec/ahc;
		e_center->re = b*h_center.re;
		e_center->im=b*h_center.im;
	 }
	return 1;
} /* h_to_e_data */

int
comp_pack_centers(p) /* find centers based on current radii.
Start with alpha vert; normalize gamma vert on y>0 axis; then use
order of faces. Only compute each circle's center once. */
struct p_data *p;
{
	int nf,i,vert;
	complex c,d;

	for (i=1;i<=p->nodecount;i++) p->packK_ptr[i].plot_flag=0;
	nf=p->first_face;
	place_face(p,nf,p->faces[nf].index_flag);
	for (i=0;i<3;i++) p->packK_ptr[p->faces[nf].vert[i]].plot_flag=1;
	while ( (nf=p->faces[nf].next_face)!=p->first_face )
	 {
		vert=p->faces[nf].vert[(p->faces[nf].index_flag +2) % 3];
		if (!p->packK_ptr[vert].plot_flag) /* not yet computed */
			comp_center_face(p,nf,-1);
		p->packK_ptr[vert].plot_flag=1;
	 }
/* normalize */
	c=p->packR_ptr[p->alpha].center;
	d=p->packR_ptr[p->gamma].center;
	if (p->hes<0) h_norm_pack(p,c,d);
	return 1;
} /* comp_pack_centers */

comp_center_face(p,f,flag) /* compute and store center/rad of last circle of
face f in pack p; flag is optional index, else use index_flag of f_data. */
int f,flag;
struct p_data *p;
{
	int j,k,v,index;
	float o1,o2,o3;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=p->packK_ptr;
	pR_ptr=p->packR_ptr;
	if (f<1 || f> p->facecount) return;
	if (flag<0 || flag>2) index=p->faces[f].index_flag;
	else index=flag;
	j=p->faces[f].vert[index];
	k=p->faces[f].vert[(index+1) % 3];
	v=p->faces[f].vert[(index+2) % 3];
	o1=o2=o3=1.0;
	h_compcenter(pR_ptr[j].center,pR_ptr[k].center,
		&pR_ptr[v].center,pR_ptr[j].rad,pR_ptr[k].rad,
		&pR_ptr[v].rad,o1,o2,o3);
} /* comp_center_face */

int
h_norm_pack(p,c,d) /* normalizes hyperbolic data of pack p by 
putting point c at origin and d on pos y-axis */
complex c,d;
struct p_data *p;
{
	int i;
	float abd;
	complex temp;
	struct R_data *pR_ptr;

	pR_ptr=p->packR_ptr;
	if (cAbs(c)>=OKERR) /* j vertex not origin */
		for (i=1;i<=p->nodecount;i++)
			pR_ptr[i].center=mob_norm(pR_ptr[i].center,c,d);
	else if (d.im>=OKERR || d.im<=(-OKERR) || d.re<=(-OKERR))
		/* just a rotation is needed */
	 {
		abd=1.0/cAbs(d);
		d.re *= abd;
		d.im *= abd;
		for (i=1;i<=p->nodecount;i++)
			pR_ptr[i].center=cdiv(pR_ptr[i].center,d);
	 }
		/* now rotate another 90 degrees --- need to clean this up
		   later */
	for (i=1;i<=p->nodecount;i++)
	 {
		temp=pR_ptr[i].center;
		pR_ptr[i].center.re=(-temp.im);
		pR_ptr[i].center.im=temp.re;
	 }
} /* h_norm_pack */

int
place_face(p,face,index) /* compute cents of face, index at origin,
next vert in standard relation */
struct p_data *p;
int face,index;
{
	int a,k;
	float r,r2,ovlp,erad,s,ss1,ss2,s1,s2;
	struct K_data *pK_ptr;
	struct R_data *pR_ptr;

	pK_ptr=p->packK_ptr;
	pR_ptr=p->packR_ptr;
	if (face>p->facecount || face<1 || index <0 || index >2) return 0;
	a=p->faces[face].vert[index];
	k=p->faces[face].vert[(index+1) % 3];
	ovlp=1.0;
	if (p->hes<0) /* hyp case */
	 {
		s1=pR_ptr[a].rad;
		s2=pR_ptr[k].rad;
		if (s1<=0)
		 {
			s1=pR_ptr[a].rad=.01;
		 }
		pR_ptr[a].center.re=0; pR_ptr[a].center.im=0;
		if (s2<=0) /* if next one is infinite radius */
		 {
			pR_ptr[k].center.re=1;pR_ptr[k].center.im=0;
			erad=(1-s1)/(1+s1);
			pR_ptr[k].rad=(-1)*(1-erad*erad)/(2.0+2.0*erad*ovlp);
		 }
		else
		 {
			ss1=s1*s1;ss2=s2*s2;
			s=exp( acosh((1/(4.0*s1*s2))*((1+ss1)*(1+ss2)
				+(1-ss1)*(1-ss2)*ovlp) ) );
			pR_ptr[k].center.re=(s-1)/(s+1);
			pR_ptr[k].center.im=0.0;
		 }
	 }
	else if (p->hes>0) /* sphere case */
	 {
			/* alpha at south pole */
		pR_ptr[a].center.re=0;pR_ptr[a].center.im=M_PI; 
			/* next out pos x-axis */
		pR_ptr[k].center.re=0.0;
		pR_ptr[k].center.im=M_PI-pR_ptr[a].rad-pR_ptr[k].rad;
	 }
	else /* eucl case */
	 {
			/* alpha at origin */
		r=pR_ptr[a].rad;
		pR_ptr[a].center.re=0;pR_ptr[a].center.im=0; 
			/* next on x-axis */
		r2=pR_ptr[k].rad;
		pR_ptr[k].center.re=sqrt(r*r+r2*r2+2*r*r2*ovlp);
		pR_ptr[k].center.im=0; 
	 }
	comp_center_face(p,face,index);
	return 1;
} /* place_face */
		
int
h_compcenter(z1,z2,z3,s1,s2,s3,o1,o2,o3) /* given 2 hyp centers/s_radii and
s_radius of third in ordered triple, return hyp center and s_radius
for third circle. oj is cos of overlap angle opposite to circle j. 
There is no consistency check on first two circles' data. */
complex z1,z2,*z3;
float s1,s2,*s3,o1,o2,o3;
{
	float flag=1,cc,rad,ahc,sc,s,ac,x0,x1,r,sidelength,ss1,ss3;
	complex a,b,c,c1,w1,w2,w3,cent2,newcent3,par;

	a=z1;
	b=z2;
	if ((s1<=0) && (s2>0)) /* second circle finite, first not */
	 {
		a=z2; b=z1; s=s1; s1=s2; s2=s; s=o1; o1=o2; o2=s; flag=(-1); 
			/* interchange order */
	 }
	if (s1>0) /* first is now finite */
	 {
		cc=h_cos_overlap(s1,s2,*s3,o1,o2,o3);
		z3->re=cc;
		z3->im=flag*sqrt(1-cc*cc);
		if (*s3>0)
		 {
			ss1=s1*s1;
			ss3=*s3*(*s3);
			sidelength=exp( acosh(
			  (1/(4.0*s1*(*s3)))*
			  ((1+ss3)*(1+ss1)+(1-ss3)*(1-ss1)*o2) ) );
			ahc=(sidelength-1)/(sidelength+1); /* abs value of hyp 
							center */
			z3->re *= ahc;
			z3->im *= ahc; /* center as if z1 at origin */
			*z3=mob_norm_inv(*z3,a,b);   /* move to right place */
			return (0);
		 }
		r=(1-s1)/(1+s1);
		sc=(r*r+1+2*r*o2)/(2*(1+r*o2));
			/* abs value of eucl center c of third circle */
		c.re=sc*z3->re;
		c.im=sc*z3->im;
		rad=1-sc; /* Now have c and its radius */
		w1.re=c.re-rad;
		w1.im=c.im;
		w2.re=c.re+rad;
		w2.im=c.im;
		w3.re=c.re;
		w3.im=c.im+rad; /* three points on the circle */
		w1=mob_norm_inv(w1,a,b);
		w2=mob_norm_inv(w2,a,b);
		w3=mob_norm_inv(w3,a,b);       /* 3 pts on new circle */
		circle_3(w1,w2,w3,&c,&rad);
	 	*s3=(-rad); /* store new s_radius */
		ac=cAbs(c);
		z3->re=c.re/ac;
		z3->im=c.im/ac; /* get hyp center on unit circle */
		return (0);
	 }
	/* remaining: first 2 have infinite rad.*/
	if (*s3<=0) /* third also infinite. As temp shortcut, we make 
		radius large but finite. */
		*s3=OKERR;
	/* Now, third finite.
		Pretend 3 is at origin, locate first/second on unit circle.
		Apply Mobius to make first correct eucl size (keeping real
		axis fixed), move back by rotation to get first to 
		original location. */
	cent2.re=h_cos_overlap(*s3,s1,s2,o3,o1,o2);
	cent2.im=sqrt(1-cent2.re*cent2.re); /* circle 2 center in 
		normalize position */
	r=(1-*s3)/(1+*s3); /* eucl radius of circle 3 when at origin */
	x0=1-(r*r+1+2*r*o2)/(2*(1+r*o2)); /* eucl center of circle 1 
		in normalized position - tangent at 1. */
	x1=1-fabs(s1); 	 /* desired eucl center of cir 1 */
	if (x1<OKERR || x1>(1-OKERR)) x1=.5; /* in case s1 
		hadn't been set */
	par.re=(x1-x0)/(x1*x0-1);
	par.im=0.0; /* set parameter for mobius 
			(z-par)/(1-conj(par)*z) */
	newcent3.re=(-1)*par.re;
	newcent3.im=0.0; /* this is just mob_trans(ORIGIN,par) */
	*z3=cmult(z1,newcent3); /* rotate to put circle 1 back in right
		location */
	return (0);	
} /* h_compcenter */

float
h_cos_overlap(s1,s2,s3,t1,t2,t3) /* given ordered triple of s_radii and
cosines of overlap angles, compute cosine of angle at the first
circle of triangle they form. Note: tj is cos of overlap angle opposite 
circle j. */
float s1,s2,s3,t1,t2,t3;
{
	float h1,h2,h3,e1,e2,e3,L,len,ans;

	if (s1<=0) return (1.0);
	e1=(1.0-s1)/(1.0+s1);
	h1=-log(s1);
	if (s2<=0) L=1.0;
	else 
	 {
		h2=-log(s2);
		len=exp(-h2-acosh(cosh(h1)*cosh(h2)+sinh(h1)*sinh(h2)*t3));
		L=(1.0-len)/(1.0+len);
	 }
	e2=(L*L-e1*e1)/(2.0*e1*t3+2.0*L);
	if (s3<=0) L=1.0;
	else 
	 {
		h3=-log(s3);
		len=exp(-h3-acosh(cosh(h1)*cosh(h3)+sinh(h1)*sinh(h3)*t2));
		L=(1.0-len)/(1.0+len);
	 }
	e3=(L*L-e1*e1)/(2.0*e1*t2+2.0*L);
		/* have euclidean radii now */
	ans=e_cos_overlap(e1,e2,e3,t1,t2,t3); /* find cos from eucl routine */
	return (ans);
} /* h_cos_overlap */

complex
mob_norm_inv(w,a,b) /* returns preimage of w under mobius of disc
which maps a to zero and b to positive x-axis */
complex w,a,b;
{
	complex z;
	float c;

	z=mob_trans(b,a);
	c=cAbs(z);
	z.re/=c;
	z.im/=c;
	z=mob_trans(cmult(w,z),a);
	return (z);
} /* mob_norm_inv */

float
e_cos_overlap(e1,e2,e3,t1,t2,t3) /* given three eucl radii and cosines
of opposite overlap angles, return cos of angle at e1 */
float e1,e2,e3,t1,t2,t3;
{
	float l1,l2,l3,ang;

	l3=e1*e1+e2*e2+2*e1*e2*t3;
	l2=e1*e1+e3*e3+2*e1*e3*t2;
	l1=e2*e2+e3*e3+2*e2*e3*t1;
	ang=(l2+l3-l1)/(2*sqrt(l2*l3));
	return (ang);
} /* e_cos_overlap */

int
circle_3(z1,z2,z3,c,rad) /* find eucl center and rad for circle
through 3 given points. Returns 1 in case of success. */
complex z1,z2,z3;
complex *c;
float *rad;
{
	float det,a1,a2,b1,b2,c1,c2,dum;

	a1=z2.re-z1.re;
	a2=z3.re-z2.re;
	b1=z2.im-z1.im;
	b2=z3.im-z2.im;
	det=2*(a1*b2-b1*a2);
	if (fabs(det)<OKERR) return 0;
	dum=cAbs(z2);
	dum=dum*dum;
	c1=cAbs(z1);
	c1=(-c1*c1);
	c1+=dum;
	c2=cAbs(z3);
	c2=c2*c2;
	c2-=dum;
	c->re=(b2*c1-b1*c2)/det;
	c->im=(a1*c2-a2*c1)/det;
	*rad=cAbs(csub(*c,z1));
	return 1;
} /* circle_3 */

/* ===================  Combinatoric processing =============== */
/* much more than needed for this limited goal */

int
complex_count(p) /* Identify and count the no. of bdry and interior
nodes, bdry and interior components,  Euler characteristic, genus, etc,
for pack p. Stores indicators of bdry and int components in packdata (at
most MAX_COMPONENTS of each). Sets default color codes for
circles/faces.*/
struct p_data *p;
{
	int i,node,count,bcount,new_vert,num_edges,m,n;
	int comp_count,ptr,k,j,fj,icount,more_flag,next_bdry;
	struct K_data *pK_ptr;

	pK_ptr=p->packK_ptr;
/* some initialization and counting */
	node=p->nodecount;
	for (i=1;i<=MAX_COMPONENTS;i++) 
		p->bdry_starts[i]=p->int_starts[i]=0;
	bcount=0;
	comp_count=0;
	for (i=1;i<=node;i++) 
	 {
		pK_ptr[i].plot_flag=0;
		if (pK_ptr[i].flower[0]!=pK_ptr[i].flower[pK_ptr[i].num]) 
		 {
			bcount++;
			pK_ptr[i].bdry_flag=1;
if (printf_bool)
printf("node %2d identified as a boundary node\n",i);
		 }
		else pK_ptr[i].bdry_flag=0;
	 }
	p->intnode=node-bcount;

/* identify bdry components, give pointers; note: index starts with 1 */
	count=0;
	new_vert=0;
	do
	 {
		new_vert++;
		if (pK_ptr[new_vert].bdry_flag && !pK_ptr[new_vert].plot_flag)
		 {
		   comp_count++;
		   if (comp_count>MAX_COMPONENTS)
		    {
			return 0;
		    }
		   p->bdry_starts[comp_count]=new_vert;
		   next_bdry=new_vert;
		   do
		    {
			pK_ptr[next_bdry].plot_flag=1;
			next_bdry=pK_ptr[next_bdry].flower[0];
			count++;
			if (count>p->nodecount)
			 {
				return 0;
			 }
		    }
		   while (next_bdry!=new_vert);
		 }
	 }
	while (count<bcount);
	p->num_bdry_comp=comp_count;

/* identify int components and pointers to them */
	icount=node-bcount; /* number of interiors */
	if (icount>0)
	 {
		comp_count=0;
		do
 {
	ptr=0;
	do {ptr++;}
	while (pK_ptr[ptr].plot_flag); 
		/* find smallest index not plotted (hence int node) */
	comp_count++;
	if (comp_count>MAX_COMPONENTS)
	 {
		return 0;
	 }
	icount--;
	p->int_starts[comp_count]=ptr;
	pK_ptr[ptr].plot_flag=ptr;
	do
	 {
	   more_flag=0; /* will tell me whether I am still finding anything */
	   for (k=ptr;k<=node;k++)
	    {
	      if (pK_ptr[k].plot_flag==ptr)
	       {
		for (j=0;j<=pK_ptr[k].num;j++) 
		 {
			fj=pK_ptr[k].flower[j];
			if (!pK_ptr[fj].plot_flag) 
			 {
				pK_ptr[fj].plot_flag=ptr;
				icount--;
				more_flag=1;
			 }
		 }
	       }
	    }
	 } /* end of inside do */
	while (more_flag && icount>0);
 } /* end of outside do */
		while (icount>0);
		p->num_int_comp=comp_count;
	 } /* end of 'if icount' */
/* find Euler characteristic and genus */
	count=0;
	for (k=1;k<=node;k++) count += pK_ptr[k].num;
	p->facecount=count/3;
	num_edges=(count + bcount)/2;
	p->euler=node-num_edges+p->facecount;
	p->genus=(2-p->euler-p->num_bdry_comp)/2;
	if (packdata.genus!=0 || packdata.euler!=1) return 0; 
		/* not simply connected */
/* allocate space for face data, store ordered triples of verts. */
	if (!alloc_faces_space(p))
	 {
		return 0;
	 }
	count=1;
	i=1;
	while (count<=p->facecount && i<=p->nodecount)
	{
		for (j=0;j<pK_ptr[i].num;j++)
		 {
			m=pK_ptr[i].flower[j];
			n=pK_ptr[i].flower[j+1];
			if (m>i && n>i) 
			 {
				p->faces[count].vert[0]=i;
				p->faces[count].vert[1]=m;
				p->faces[count].vert[2]=n;
				count++;
			 }
		 }
		i++;
	 }
	return 1;
} /* complex_count */

int
alloc_faces_space(p) /* allocate face data space (knowing facecount) */
struct p_data *p;
{
	if (p->faces!=NULL) free(p->faces);
	p->faces=(f_data *)NULL;
	p->faces=(f_data *)malloc((p->facecount+1)*sizeof(f_data));
	if ((p->faces)==NULL) return 0;
	return 1;
} /* alloc_faces_space */

int
facedraworder(p,drawlist) /* return 1 if succeed. drawlist not yet used */
struct p_data *p;
struct Vertlist *drawlist;
{
	int i,j,nfaces,loop_count=0,count=1,alpha,*face_ptr,f,num;
	int tripflag,stopface,dummy;
	int stopvert,lastvert,stop=0,looplimit,stopflag;
	struct stface *fhandle,*chandle,*temp,*killit;
	struct Vertlist *face_order,*fdo;

/* initialize */

	alpha=p->alpha;
	build_org(p);
	for (i=1;i<=p->facecount;i++)
	 {
		p->faces[i].rwb_flag=2;
		p->faces[i].plot_flag=0;
	 }
	chandle=temp=(struct stface *)malloc(sizeof(struct stface));
	face_order=fdo=(struct Vertlist *)malloc(sizeof(struct Vertlist));

/* start with faces about alpha */

	face_ptr=face_org+(alpha*(NUM_PETALS+1))+1;
	num=p->packK_ptr[alpha].num;
	chandle->face=f=*(face_ptr); /* first face */
	fdo->v=f;
	p->faces[f].plot_flag=1;
	p->faces[f].rwb_flag=1;
	p->faces[f].index_flag=find_index(p,f,alpha);
	for (i=1;i<num;i++)
	 {
		f=*(face_ptr+i);
		chandle->next=(struct stface *)malloc(sizeof(struct stface));
		fhandle=chandle;
		chandle=chandle->next;
		chandle->prev=fhandle;
		chandle->face=f;
		fdo=fdo->next=
			(struct Vertlist *)malloc(sizeof(struct Vertlist));
		fdo->v=f;
		p->faces[f].plot_flag=1;
		p->faces[f].rwb_flag=1;
		p->faces[f].index_flag=find_index(p,f,alpha);
	 }
	fdo->next=NULL;
	chandle->next=temp;
	temp->prev=chandle;
	fdo->next=NULL;

	/* -------------- the rest (main loop) -------------------- */
	looplimit=(100)*p->nodecount;
	while (!stop && count && loop_count<looplimit)
	 {
		count=stopface=0;
		while ( (stopflag=
		   zip_vert(p,&chandle,&fdo,stopface,&stopface,&tripflag))>0 
		   && loop_count<looplimit)
			/* go until you have a marker vert */
		 {loop_count++;count += stopflag;}
		if (stopflag<0) stop=1;/* redchain a closed flower */
		while ( !stop && (stopflag=
		   zip_vert(p,&chandle,&fdo,stopface,&dummy,&tripflag))>=0 
		   && !tripflag && loop_count<looplimit)
			/* go until you encounter stopvert or stop */
		 {loop_count++;count += stopflag;}
		if (stopflag<0) stop=1; /* redchain a closed flower */
		while ( !stop && (stopflag=
		   zip_vert(p,&chandle,&fdo,stopface,&stopface,&tripflag))>0  
		   && loop_count<looplimit)
			/* go until you have a marker vert */
		 {loop_count++;count += stopflag;}
		if (stopflag<0 || stop) /* redchain a closed flower */
		 {
			stop=1;
			p->faces[chandle->face].rwb_flag=0;
			killit=chandle=chandle->next;
			while (chandle->next!=killit)
			 {
				p->faces[chandle->next->face].rwb_flag=0;
				chandle=chandle->next;
			 }
			free_face(&chandle);
			chandle=NULL;
			p->first_red_face=0;
		 }
	 } /* end of while */
/* find blue faces to be put in red chain */
	if (!stop) add_blue(p,&chandle,&fdo,&stopvert);
	lastvert=0;
	while (!stop && add_blue(p,&chandle,&fdo,&lastvert)>=0
		&& stopvert!=lastvert );
/* arrange final drawing order */
	fdo=face_order;
	if ( (f=final_order(p,&fdo,chandle)) ) 
		/* have some missing faces */
	 {
		if (f==-1) /* error, redo all by old-style drawing order */
		 {
			for (i=1;i<=p->nodecount;i++) 
				p->packK_ptr[i].plot_flag=0;
			for (i=1;i<=p->facecount;i++) 
				p->faces[i].plot_flag=0;
			p->packK_ptr[alpha].plot_flag=1;
			f=*(face_org+(alpha*(NUM_PETALS+1))+1);
			for (j=0;j<3;j++) 
			   p->packK_ptr[p->faces[f].vert[j]].plot_flag=1;
		 } 
		if (!wrapup_order(p,fdo->v))
		 {
			return 0;
		 }
	 }
/* arrange list of red faces */	
	if (chandle==NULL) /* e.g., sphere; throw away red chain */
	 {
		for (i=1;i<=p->facecount;i++) p->faces[i].rwb_flag=0;
		p->first_red_face=0;
	 }
	else
	 {
		for (i=1;i<=p->facecount;i++) 
			p->faces[i].next_red=chandle->face;
		fhandle=chandle;
		while (chandle->next!=fhandle)
		 {
			p->faces[chandle->face].next_red=chandle->next->face;
			chandle=chandle->next;
		 }
	 }
	free_face(&chandle);
	vert_free(&face_order);
	free(face_org);
	return 1;
} /* facedraworder */

int
zip_vert(p,handle,fdo,stopface,lastface,tripflag) /* find vert and 
fan of contig red faces, reset handle downstream. If rest of faces 
blue, then zip: adjust redchain, face colors, drawing order, etc. 
Return 1 if zipped. Return 0 if not zip-able. Return -1
if redchain is around single vert. set tripflag if stopface encountered. */
struct p_data *p;
struct stface **handle;
struct Vertlist **fdo;
int stopface,*lastface,*tripflag;
{
	int cface,n,f,findex=0,bindex,eindex,i,j,num,*face_ptr,vert;
	struct stface *bzip,*new_seg,*tmp,*hold,*killit;

	bzip=*handle;
	*tripflag=0;
	cface=(*handle)->face;
	if ( (n=nghb_tri(p,(*handle)->next->face,cface))<0 ) 
	 {*lastface=cface;return 0;} 
		/* this would imply error in red chain */
	vert=p->faces[cface].vert[n];
	num=p->packK_ptr[vert].num;
	face_ptr=face_org+((vert)*(NUM_PETALS+1))+1;
/* find index of cface */
	while (*(face_ptr+findex)!=cface && findex<(num-1)) findex++;
/* find index of downstream face */
	i=findex;
	while ( (*handle)->next->face==*(face_ptr+((i+num-1) % num))
	   && (*handle)->next->face!=cface ) 
	 {
		*handle=(*handle)->next;
		i--;
	 }
	if ((*handle)->next->face==cface) return -1; 
		/* redchain goes around single vert */
	eindex= (num+i) % num;
/* find index of upstream face */
	i=findex;
	while ( bzip->prev->face==*(face_ptr+((i+1) % num))
	   && bzip->prev->face!=cface )
	 {
		bzip=bzip->prev;
		i++;
	 }
	if (bzip->prev->face==cface) return -1; 
		/* this should have been caught before */
	bindex=(i % num);
/* did we come across stopface? */
	n=(bindex+num-eindex) % num;
	for (j=0;j<=n;j++) 
		if (*(face_ptr+(j+eindex) % num)==stopface) *tripflag=1;
	*lastface=*(face_ptr+bindex);
/* only zip interiors */
	if (p->packK_ptr[vert].bdry_flag) return 0;
	i=bindex+1;
	while ( (j=(i % num))!=eindex )
	 {
		if (p->faces[*(face_ptr+j)].rwb_flag!=2) return 0;
		i++;
	 }

/* success: push redchain over vertex */

	if ( (i=(bindex+1) % num)!=eindex ) /* new red segment */
	 {
		tmp=new_seg=(struct stface *)malloc(sizeof(struct stface));
		new_seg->prev=bzip;
		new_seg->face=*(face_ptr+i);
		p->faces[new_seg->face].rwb_flag=1;
		add_face_order(p,fdo,new_seg->face,bzip->face);
		i++;
		while ( (j=i % num)!=eindex)
		 {
			hold=tmp;
			tmp=tmp->next=
			   (struct stface *)malloc(sizeof(struct stface));
			tmp->prev=hold;
			tmp->face=*(face_ptr+j);
			p->faces[tmp->face].rwb_flag=1;
			add_face_order(p,fdo,tmp->face,hold->face);
			i++;
		 }
	 }
	else /* no new red */
	 {
		tmp=bzip;
		new_seg=*handle;
	 }
/* remove old segment (becomes white) */
	killit=bzip->next;
	while (killit!=*handle)
	 {
		hold=killit->next;
		p->faces[killit->face].rwb_flag=0;
		free(killit);
		killit=hold;
	 }		
/* insert new red segment */
	(*handle)->prev=tmp;
	tmp->next=(*handle);
	bzip->next=new_seg;
	return 1;
} /* zip_vert */

int
add_blue(p,handle,fdo,vert) /* see if blue face can be added to face
in red chain. Return -1 when we've gone around until we hit a spot 
where blue face has been added (this should mean we have added all
possible blue faces. Return 1 if face is added, zero otherwise. */
struct p_data *p;
struct stface **handle;
struct Vertlist **fdo;
int *vert;
{
	int cface,n,f,findex=0,nindex,i,j,num,newface,*face_ptr;
	struct stface *hold;

	cface=(*handle)->face;
	if ( (n=nghb_tri(p,(*handle)->next->face,cface))<0 ) return -1; 
		/* this is recently added face; means we've gone all the
			way around  */
	*vert=p->faces[cface].vert[n];
	num=p->packK_ptr[*vert].num;
	face_ptr=face_org+((*vert)*(NUM_PETALS+1))+1;
	while (*(face_ptr+findex)!=cface && findex<(num-1)) findex++;
	i=findex;
	while ( (*handle)->next->face==*(face_ptr+((i+num-1) % num))
	   && (*handle)->next->face!=cface ) 
	 {
		*handle=(*handle)->next;
		i--;
	 }
	nindex=(num+i-1) % num;
	newface=*(face_ptr+nindex);
	if ( p->faces[newface].rwb_flag!=2
	   || (p->packK_ptr[*vert].bdry_flag && nindex>findex) ) 
		return 0; /* face not blue or on wrong end of bdry flower */
	hold=(*handle)->next;
	(*handle)->next=(struct stface *)malloc(sizeof(struct stface));
	(*handle)->next->prev=*handle;
	*handle=(*handle)->next;
	(*handle)->face=newface;
	(*handle)->next=hold;
	hold->prev=*handle;
	p->faces[newface].rwb_flag=1;
	add_face_order(p,fdo,newface,(*handle)->prev->face);
	(*handle)=(*handle)->next;
	return 1;
} /* add_blue */

add_face_order(p,fdo,face,preface)  /* Insert face in draw order next 
to preface. */
struct p_data *p;
struct Vertlist **fdo;
int face,preface;
{
	int n;

	(*fdo)=(*fdo)->next=(struct Vertlist*)malloc(sizeof(struct Vertlist));
	(*fdo)->v=face;
	(*fdo)->next=NULL;
	p->faces[face].plot_flag=1;
	n=nghb_tri(p,preface,face);
	if (n<0) p->faces[face].index_flag=0; /* shouldn't happen */ 
	else p->faces[face].index_flag=n;
} /* add_face_order */

int 
final_order(p,fdo,handle)  /* strategy: Simply connected case, 
do in order without regard to face color. Non-simply connected case,
start with alpha, proceed to first red face, do whole red chain,  
come back for rest of interior. This is needed 
because some circles will be moved while drawing red faces.
Return 0 for success and all faces plotted, -1 if an error occurs
which will require use of old method, 1 if success, but some faces 
still not plotted. */
struct p_data *p;
struct Vertlist **fdo;
struct stface *handle;
{
	int count=0,i,j,blue_count,f,start_face,index,fhold;
	struct Vertlist *start,*temp;
	struct stface *hold;

	start=temp=*fdo;

/* set next_face */

	p->first_face=(*fdo)->v;
	for (i=1;i<=p->nodecount;i++) p->packK_ptr[i].plot_flag=0;
	for (i=1;i<=p->facecount;i++) p->faces[i].plot_flag=0;
	for (j=0;j<3;j++) 
		p->packK_ptr[p->faces[(*fdo)->v].vert[j]].plot_flag=1;
	p->faces[p->first_face].plot_flag=1;

/* ----------------- simply connected ------------------------- */

	if ((p->euler==1 || p->euler==2) && p->genus==0)
	 {
		f=(*fdo)->v;
		p->faces[f].next_face=(*fdo)->next->v;
		*fdo=(*fdo)->next; /* take care of first face first */
		while ((*fdo)->next!=NULL) /* get all rest in fdo */
		 {
			f=(*fdo)->v;
			index=p->faces[f].index_flag;
			p->faces[f].next_face=(*fdo)->next->v;
			p->faces[f].index_flag=nice_index(p,f,index);
			p->packK_ptr[
				p->faces[f].vert[(p->faces[f].index_flag+2) % 3]
				].plot_flag=1;
			p->faces[f].plot_flag=1;
			*fdo=(*fdo)->next;
		 }
/* set first_red_face */
		if (handle!=NULL)
		 {
			while (temp->next!=NULL 
			   && p->faces[temp->v].rwb_flag!=1)
				temp=temp->next;
			if (temp->next!=NULL) p->first_red_face=temp->v;
			else p->first_red_face=0;
		 }
/* check for remaining faces */
		count=0;
		for (i=1;i<=p->facecount;i++) 
			if (p->faces[i].plot_flag) count++;
		if (count<p->facecount) return 1;
		p->faces[(*fdo)->v].next_face=start->v;
		return 0;
	 }

/* ------------------------ non-simply connected ----------------- */

	else
	 {
/* first interiors, from alpha */
		while (!p->faces[(*fdo)->v].rwb_flag && (*fdo)->next!=NULL)
		 {
			f=(*fdo)->v;
			index=p->faces[f].index_flag;
			p->faces[f].next_face=(*fdo)->next->v;
			p->faces[f].index_flag=nice_index(p,f,index);
			p->packK_ptr[
			   p->faces[f].vert[(p->faces[f].index_flag+2) % 3]
			   ].plot_flag=1;
			p->faces[f].plot_flag=1;
			*fdo=(*fdo)->next;
		 }
		if ((*fdo)->next==NULL) /* should not happen */
		 {
			p->faces[(*fdo)->v].next_face=start->v;
			return -1;
		 }
/* now, put whole red chain in draw order */
		p->first_red_face=f=(*fdo)->v;
		in_handle(&handle,(*fdo)->v);
		blue_count=0;
		while (handle->next->face!=(*fdo)->v 
			&& blue_count<=p->facecount)
		 {
			p->faces[handle->face].next_face=handle->next->face;
			handle=handle->next;
			blue_count++;
		 }
		if (blue_count>p->facecount)
		 {
			return -1;
		 }
/* loop thru again to reset index_flags. */
		in_handle(&handle,(*fdo)->v);
		start_face=f=handle->face;
		index=p->faces[f].index_flag;
		p->faces[f].index_flag=nice_index(p,f,index);
		p->packK_ptr[p->faces[f].vert[(p->faces[f].index_flag+2) % 3]
			].plot_flag=1;
		p->faces[f].plot_flag=1;
		while(handle->next->face!=start_face)
		 {
			blue_count=0;
			hold=handle;
			if ((index=nghb_tri(p,hold->face,
			   handle->next->face))<0 
			   && (index=nghb_tri(p,hold->prev->face,
			   handle->next->face))<0) return -1;
				/* blue face, so back up; if previous
				face isn't contiguous, error. */ 
			f=handle->next->face;
			p->faces[f].index_flag=nice_index(p,f,index); 
			p->packK_ptr[p->faces[f].
				vert[(p->faces[f].index_flag+2) % 3]
				].plot_flag=1;
			p->faces[f].plot_flag=1;
			handle=handle->next;
		 } /* end of while */
/* get rest of white in fdo */
		fhold=handle->face;
		while ((*fdo)->next!=NULL)	
		 {
			while ( p->faces[(*fdo)->v].rwb_flag 
			   && (*fdo)->next!=NULL )
				*fdo=(*fdo)->next;
					/* pass up any red/blue faces */
			if ((*fdo)->next!=NULL)
			 {
				p->faces[fhold].next_face=f=(*fdo)->v;
				p->faces[f].index_flag=
					nice_index(p,f,p->faces[f].index_flag);
				p->packK_ptr[p->faces[f].
				   vert[(p->faces[f].index_flag+2) % 3]
				   ].plot_flag=1;
				p->faces[f].plot_flag=1;
				fhold=(*fdo)->v;
				*fdo=(*fdo)->next;
			 }
		 }
/* check for remaining faces */
		count=0;
		for (i=1;i<=p->facecount;i++) 
			if (p->faces[i].plot_flag) count++;
		if (count<p->facecount) return 1;
		p->faces[fhold].next_face=start->v;
		return 0;
	 } /* end of non-simply connected case */
 } /* final_order */

int
wrapup_order(p,lastface) /* Wrap up face order. Put in remaining faces.
return 0 if fails. */
struct p_data *p;
int lastface;
{
	int k=1,hit=1,stop=0,F,i,index,count=0,n;

	F=p->facecount;
	for (i=1;i<=F;i++) if (p->faces[i].plot_flag) count++;
	while (count<F && !stop)
	 {
		if (k > F)
		 {
			k=1;
			if (!hit) stop=1; /* no new faces being added */
			hit=0;
		 }
		while (p->faces[k].plot_flag && k < F) k++;
		if (!p->faces[k].plot_flag && face_ready(p,k,&index))
		 {
			if ( (n=nghb_tri(p,lastface,k))<0 )
				p->faces[k].index_flag=index;
			else p->faces[k].index_flag=n;
			p->faces[lastface].next_face=k;
			p->faces[k].plot_flag=1;
			p->packK_ptr[p->faces[k].
				vert[(index+2) % 3]].plot_flag=1;
			lastface=k;
			count++;
			hit=1;
		 }
		k++;
	 }
	p->faces[lastface].next_face=p->first_face;
	if (count<F) return 0; /* shouldn't ever happen */
	return 1;
} /* wrapup_order */

int 
face_ready(p,face,index) /* return 1 if face ready to draw, set
index of vert to draw first. */
struct p_data *p;
int face,*index;
{
	if (p->packK_ptr[p->faces[face].vert[0]].plot_flag 
		&& p->packK_ptr[p->faces[face].vert[1]].plot_flag) 
	 {
		*index=0;
		return 1;
	 }
	if (p->packK_ptr[p->faces[face].vert[1]].plot_flag 
		&& p->packK_ptr[p->faces[face].vert[2]].plot_flag) 
	 {
		*index=1;
		return 1;	
	 }
	if (p->packK_ptr[p->faces[face].vert[2]].plot_flag 
		&& p->packK_ptr[p->faces[face].vert[0]].plot_flag) 
	 {
		*index=2;
		return 1;
	 }
	*index=0;
	return 0;
} /* face_ready */

in_handle(handle,face) /* position closed list at face */
struct stface **handle;
int face;
{
	struct stface *current;

	current=*handle;
	if ((*handle)->face==face) return;
	do {current=current->next;}
	while ( (current->face!=face) && ((current->next)!=*handle) );
	*handle=current;
} /* in_handle */

int
nice_index(p,f,index) /* set index_flag of face f so that it uses nghb'ing
white face using verts already drawn, if possible. "index" is first 
preference and default. plot_flags are set for verts already drawn. */
struct p_data *p;
int f,index;
{
	int v,w,i,m,num,g,*face_ptr;

	for (i=0;i<3;i++)
	 {
		v=p->faces[f].vert[(index+i)%3];
		w=p->faces[f].vert[(index+i+1)%3];
		if (p->packK_ptr[v].plot_flag && !p->packK_ptr[v].bdry_flag
		   && p->packK_ptr[w].plot_flag && !p->packK_ptr[w].bdry_flag)
		 {
			face_ptr=face_org+(w*(NUM_PETALS+1));
			num=*face_ptr;
			face_ptr++;
			m=0;
			while (*(face_ptr+m)!=f && m<(num-1)) m++;
			g=*(face_ptr+((m+1) % num));
			if (p->faces[g].rwb_flag==0) return ((index+i)%3);
				 /* yes, face is white */
		 }
	 }
	for (i=0;i<3;i++) /* fallback: try for white neighbor */
	 {
		v=p->faces[f].vert[(index+i)%3];
		w=p->faces[f].vert[(index+i+1)%3];
		if (p->packK_ptr[v].plot_flag 
		   && p->packK_ptr[w].plot_flag)
		 {
			face_ptr=face_org+(w*(NUM_PETALS+1));
			num=*face_ptr;
			face_ptr++;
			m=0;
			while (*(face_ptr+m)!=f && m<(num-1)) m++;
			g=*(face_ptr+((m+1) % num));
			if (p->faces[g].rwb_flag==0) return ((index+i)%3);
				 /* yes, face is white */
		 }
	 }
	for (i=0;i<3;i++) /* fallback: use any two plotted ngbh's. */
	 {
		v=p->faces[f].vert[(index+i)%3];
		w=p->faces[f].vert[(index+i+1)%3];
		if (p->packK_ptr[v].plot_flag 
		   && p->packK_ptr[w].plot_flag) return ((index+i)%3);
	 }
	if (index<0 || index>2) return 0;
	return index;
} /* nice_index */

free_face(handle)
struct stface **handle;
{
	struct stface *kill;

	if ((*handle)==NULL) return;
	while ((*handle)->next!=*handle) 
	 {
		(*handle)->prev->next=(*handle)->next;
		(*handle)->next->prev=(*handle)->prev;
		kill=*handle;
		*handle=(*handle)->prev;
		free(kill);
	 }
	free(*handle);
} /* free_face */

int
build_org(p) /* create array, row for each vertex list which faces it's in */
struct p_data *p;
{
	int i,j,vert,nodecount,*number,num,v,w,spot,*face_ptr,f,hold;

	nodecount=p->nodecount;
	if ( (face_org=(int *)malloc((NUM_PETALS+1)*
		(nodecount+1)*sizeof(int)) )==NULL )
	 {
		return 0;
	 }
	for (j=1;j<=nodecount;j++) *(face_org+(j*(NUM_PETALS+1)))=0; 
	for (i=1;i<=p->facecount;i++) 
	for (j=0;j<3;j++)
	 {
		vert=p->faces[i].vert[j];
		number=face_org+(vert*(NUM_PETALS+1));
		*(number+(*number)+1)=i;
		*(number) +=1;
	 }
/* put in order */
	for (j=1;j<=nodecount;j++)
	 {
		face_ptr=face_org+(j*(NUM_PETALS+1));
		num=p->packK_ptr[j].num;
		for (i=1;i<=num;i++)
		 {
			w=p->packK_ptr[j].flower[i-1];
			spot=i;
			while (spot<=num 
			   && (f=check_face(p,*(face_ptr+spot),j,w))==(-1) ) 
				spot++;
			hold=*(face_ptr+i);
			*(face_ptr+i)=*(face_ptr+spot);
			*(face_ptr+spot)=hold;
		 }
	 }
	return 1;				
} /* build_org */

/* ======================= utility routines ================ */

complex
mob_norm(z,a,b) /* returns value at z of mobius transform of unit
disc which maps a to zero and b to positive x-axis. */
complex z,a,b;
{
	complex w,v;
	float c;

	v=mob_trans(b,a);
	c=cAbs(v);
	if (c<OKERR) return (mob_trans(z,a));
	w=cdiv(mob_trans(z,a),v);
	w.re*=c;
	w.im*=c;
	return (w);
} /* mob_norm */

complex
mob_trans(z,a) /* return value for (a-z)/(1-z*conj(a)) */
complex z,a;
{
	complex w,COMP_UNIT;

	COMP_UNIT.re=1.0;COMP_UNIT.im=0.0;
	w=cdiv(csub(a,z),csub(COMP_UNIT,cmult(z,cconj(a)) ) );
	return (w);
} /* mob_trans */

complex cadd(z1,z2) complex z1,z2;
{
	complex z;
	z.re=z1.re+z2.re;
	z.im=z1.im+z2.im;
	return (z);
}

complex cmult(z1,z2) complex z1,z2;
{
	complex z;
	z.re=z1.re*z2.re-z1.im*z2.im;
	z.im=z1.re*z2.im+z2.re*z1.im;
	return (z);
}

complex csub(z1,z2) complex z1,z2;
{
	complex z;
	z.re=z1.re-z2.re;
	z.im=z1.im-z2.im;
	return (z);
}

complex cdiv(z1,z2) complex z1,z2;
{
	complex z;
	float av;

	av = 1/(cAbs(z2));
	z=cmult(z1,cconj(z2));
	z.re=av*av*z.re;
	z.im=av*av*z.im;
	return (z);
}

complex cconj(z1) complex z1;
{
	complex z;
	z.re=z1.re;
	z.im=(-z1.im);
	return (z);
}

float cAbs(z) complex z;
{
	float a;
	a=(float)sqrt((z.re)*(z.re)+(z.im)*(z.im));
	return (a);
} 

vert_free(vertlist)
struct Vertlist **vertlist;
{
	struct Vertlist *trace,*clobber;

	if (*vertlist==NULL) return;
	trace=*vertlist;
	while (trace!=NULL)
	 {
		clobber=trace;
		trace=trace->next;
		free(clobber);
	 }
	*vertlist=NULL;
} /* vert_free */

int
nghb_tri(p,m,n) /* -1 if face n doesn't share edge e with m; else gives 
index of begin vertex of e in data of n. */
int m,n;
struct p_data *p;
{
	int nj,mj,v1,v2;
	struct K_data *pK_ptr;

	pK_ptr=p->packK_ptr;
	if (m<1 || m > p->facecount || n<1 || n > p->facecount) return -1;
	for (nj=0;nj<=2;nj++)
	 {
		v1=p->faces[n].vert[nj];
		v2=p->faces[n].vert[(nj+1)%3];
		for (mj=0;mj<=2;mj++)
			if ( (v1==p->faces[m].vert[mj]) && 
			 (v2==p->faces[m].vert[(mj+2)%3]) )
				return nj;
	 }
	return -1;
} /* nghb_tri */


int
check_face(p,f,v,w) /* Return -1 if ordered edge <v,w> not in face f. 
Else, return index of v in verts of f. */
struct p_data *p;
int f,v,w;
{
	int m=0,ans=0;

	while (m<3)
	 {
		if (p->faces[f].vert[m]==v && p->faces[f].vert[(m+1)%3]==w)
			return m; 
		m++;
	 }
	return -1;
} /* check_face */
	
int
find_index(p,f,v) /* return -1 if v not in face f; else return index of v
in vertices of f. */
struct p_data *p;
int f,v;
{
	int m=0;

	while (m<3)
	 {
		if (p->faces[f].vert[m]==v) return m;
		m++;
	 }
	return -1;
} /* find_index */


