#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 1000
#define BINS 20

#define RND drand48()
#define PI 3.1416

// from https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
double pDistance(double x, double y, double x1, double y1, double x2, double y2) {

  double A = x - x1;
  double B = y - y1;
  double C = x2 - x1;
  double D = y2 - y1;

  double dot = A * C + B * D;
  double len_sq = C * C + D * D;
  double param = -1;
  if (len_sq != 0) //in case of 0 length line
    param = dot / len_sq;

  double xx, yy;

  if (param < 0) {
    xx = x1;
    yy = y1;
  }
  else if (param > 1) {
    xx = x2;
    yy = y2;
  }
  else {
    xx = x1 + param * C;
    yy = y1 + param * D;
  }

  double dx = x - xx;
  double dy = y - yy;
  return sqrt(dx * dx + dy * dy);
}

// from https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
char get_line_intersection(float p0_x, float p0_y, float p1_x, float p1_y, 
			   float p2_x, float p2_y, float p3_x, float p3_y, float *i_x, float *i_y)
{
  float s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

  float s, t;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
      // Collision detected
      if (i_x != NULL)
	*i_x = p0_x + (t * s1_x);
      if (i_y != NULL)
	*i_y = p0_y + (t * s1_y);
      return 1;
    }

  return 0; // No collision
}

int main(void)
{
  double *IR, *PA;
  double theta, theta2;
  double tx1, ty1, tx2, ty2;
  int i, j;
  int cx, cy;
  double x, y;
  double x0[NMAX], x1[NMAX], x2[NMAX], y0[NMAX], y1[NMAX], y2[NMAX];
  int nbranch;
  double r = 4;
  int DIM = 150;
  int NPLASTID = 150;
  int MAXL = 6;
  int collision, thiscollision;
  double countcollision;
  double d;
  double thisd;
  double RHO = 10;
  float ix, iy;
  FILE *fp;
  char fstr[100];
  int expt;
  double L = 10;
  double histIR[BINS+1], histPA[BINS+1];
  double dtheta = 0.1, dtheta2 = 0.5;
  double dist;
  int sample, nsample;
  int idist;
  int plastid;
  
  IR = (double*)malloc(sizeof(double)*DIM*DIM);
  PA = (double*)malloc(sizeof(double)*DIM*DIM);

  srand48(12);

  expt = 1;
  for(i = 0; i < BINS+1; i++)
    histIR[i] = histPA[i] = 0;
  nbranch = 0;
	  
  for(plastid = 0; plastid < NPLASTID; plastid++)
    {
      x0[plastid] = RND*DIM-DIM/2; y0[plastid] = RND*DIM-DIM/2;
	      
      switch(expt) {
	// no stromules
      case 0:  break;
	// single stromules
      case 1: 
	theta = RND*3.14159*2; x0[plastid] = RND*DIM-DIM/2; y0[plastid] = RND*DIM-DIM/2; r = (x0[plastid] < 0 ? 0 : 4); L = (x0[plastid] < 0 ? 0.01 : RND*MAXL);
	x1[nbranch+0] = x0[plastid]; y1[nbranch+0] = y0[plastid]; x2[nbranch+0] = x0[plastid]+(r+L)*cos(theta); y2[nbranch+0] = y0[plastid]+(r+L)*sin(theta);
	nbranch++;
	break;
	// branched stromules
      case 2: 
	theta = RND*3.14159*2; x0[plastid] = RND*DIM-DIM/2; y0[plastid] = RND*DIM-DIM/2;
	x1[nbranch+0] = x0[plastid]; y1[nbranch+0] = y0[plastid]; x2[nbranch+0] = x0[plastid]+(r+L/3)*cos(theta); y2[nbranch+0] = y0[plastid]+(r+L/3)*sin(theta);
	x1[nbranch+1] = x2[nbranch+0]; y1[nbranch+1] = y2[nbranch+0]; x2[nbranch+1] = x1[nbranch+1]+(L/3)*cos(theta+PI+2*PI/3); y2[nbranch+1] = y1[nbranch+1]+(L/3)*sin(theta+PI+2*PI/3);
	x1[nbranch+2] = x2[nbranch+0]; y1[nbranch+2] = y2[nbranch+0]; x2[nbranch+2] = x1[nbranch+2]+(L/3)*cos(theta+PI-2*PI/3); y2[nbranch+2] = y1[nbranch+2]+(L/3)*sin(theta+PI-2*PI/3);
	nbranch += 3;
	break;
	// double stromules
      case 3: 
	theta = RND*3.14159*2; x0[plastid] = RND*DIM-DIM/2; y0[plastid] = RND*DIM-DIM/2;
	x1[nbranch+0] = x0[plastid]; y1[nbranch+0] = y0[plastid]; x2[nbranch+0] = x0[plastid]+(r+L/2)*cos(theta); y2[nbranch+0] = y0[plastid]+(r+L/2)*sin(theta);
	x1[nbranch+1] = x0[plastid]; y1[nbranch+1] = y0[plastid]; x2[nbranch+1] = x0[plastid]+(r+L/2)*cos(theta+PI); y2[nbranch+1] = y0[plastid]+(r+L/2)*sin(theta+PI);
	nbranch += 2;
	break;
      }
    }
  r = 1;
  fp = fopen("geom-snaps-vmany.csv", "w");
  fprintf(fp, "expt,dist,sample,x1,y1,x2,y2\n");
  for(i = 0; i < nbranch; i++)
    {
      fprintf(fp, "%i,%f,%i,%f,%f,%f,%f\n", expt, dist, sample, x1[i], y1[i], x2[i], y2[i]);
    }
  fclose(fp);
	
  fp = fopen("vmany-grid.csv", "w");
  fprintf(fp, "x,y,ir,pa\n");
	
  for(cx = 0; cx < DIM; cx++)
    {
      printf("%i\n", cx);
      for(cy = 0; cy < DIM; cy++)
	{
	  x = cx-DIM/2; y = cy-DIM/2;
	  d = DIM*DIM;
	  for(i = 0; i < NPLASTID; i++)
	    {
	      thisd = sqrt((x-x0[i])*(x-x0[i]) + (y-y0[i])*(y-y0[i])) - r;
	      if(thisd < d) d = thisd;
	    }
	  // check distances to any stromule segments to see if these are shorter
	  for(i = 0; i < nbranch; i++)
	    {
	      thisd = pDistance(x, y, x1[i], y1[i], x2[i], y2[i]);
	      if(thisd < d || i == 0) d = thisd;
	    }
	  IR[cy*DIM+cx] = d;

	  countcollision = 0;
	  for(theta = 0; theta < 2*3.1415; theta += dtheta)
	    {
	      thiscollision = 0;
	      for(i = 0; i < nbranch && thiscollision == 0; i++)
		{
		  // check for collision between this stromule branch and a line extended at angle theta from the current point outside the boundary of the cell
		  collision = get_line_intersection(x, y, x+2*DIM*sin(theta), y+2*DIM*cos(theta), x1[i], y1[i], x2[i], y2[i], &ix, &iy);
		  if(collision)
		    {
		      // if the collision point lies within distance RHO of this point, record the collision
		      if((x-ix)*(x-ix) + (y-iy)*(y-iy) < RHO*RHO)
			{
			  thiscollision = 1;
			  break;
			}
		    }
		}

	      // build the plastid circumference out of many small line segments, and look for collisions with any of these
	      for(theta2 = 0; theta2 < 2*3.1415 && thiscollision == 0; theta2 += dtheta2)
		{
		  for(j = 0; j < NPLASTID; j++)
		    {
		      // small segment for plastid j
		      tx1 = x0[j] + r*sin(theta2); ty1 = y0[j] + r*cos(theta2);
		      tx2 = tx1 + dtheta2*r*sin(theta2-PI/2); ty2 = ty1 + dtheta2*r*cos(theta2-PI/2);
		      // detect collision
		      collision = get_line_intersection(x, y, x+2*DIM*sin(theta), y+2*DIM*cos(theta), tx1, ty1, tx2, ty2, &ix, &iy);
		      if(collision)
			{
			  // if collision point lies within RHO, record it
			  if((x-ix)*(x-ix) + (y-iy)*(y-iy) < RHO*RHO)
			    {
			      thiscollision = 1;
			      //printf("plastid collision\n");
			      break;
			    }
			}
		    }

		}
			      
	      countcollision += thiscollision;
	    }
	  PA[cy*DIM+cx] = (double)countcollision/(theta/dtheta);
	  fprintf(fp, "%i,%i,%f,%f\n", cx-DIM/2, cy-DIM/2, IR[cy*DIM+cx], PA[cy*DIM+cx]);
	}
    }
  fclose(fp);		

	    
		
  return 0;
}
  
		  
