#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NMAX 10
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
  int i;
  int cx, cy;
  double x, y;
  double x1[NMAX], x2[NMAX], y1[NMAX], y2[NMAX];
  int NBRANCH;
  double r = 4;
  int DIM = 200;
  int collision, thiscollision;
  double countcollision;
  double d;
  double thisd;
  double RHO = 10;
  float ix, iy;
  FILE *fp;
  char fstr[100];
  int expt;
  double L = 20;
  double histIR[BINS+1], histPA[BINS+1];
  double dtheta = 0.1, dtheta2 = 0.2;
  double dist;
  int sample, nsample;
  double x10, y10, x20, y20;
  int idist;
  
  IR = (double*)malloc(sizeof(double)*DIM*DIM);
  PA = (double*)malloc(sizeof(double)*DIM*DIM);

  fp = fopen("hists-many.csv", "w");
  fprintf(fp, "expt,dist,hist,val,n\n");
  fclose(fp);
  fp = fopen("geom-snaps.csv", "w");
  fprintf(fp, "expt,dist,sample,x1,y1,x2,y2\n");
  fclose(fp);
  
  for(expt = 0; expt <= 3; expt++)
    {
      for(idist = 0; idist <= 5; idist++)
	{
	  switch(idist) {
	  case 0: dist = 0; break;
	  case 1: dist = 10; break;
	  case 2: dist = 20; break;
	  case 3: dist = 30; break;
	  case 4: dist = 40; break; 
	  case 5: dist = 50; break;
	  }
	  
	  for(i = 0; i < BINS+1; i++)
	    histIR[i] = histPA[i] = 0;

	  if(expt == 0) nsample = 1;
	  else nsample = 10;
	  
	  for(sample = 0; sample < nsample; sample++)
	    {
	      printf("%i %f %i\n", expt, dist, sample);
	      x10 = 0; y10 = 0;
	      x20 = dist; y20 = 0;
	      switch(expt) {
		// no stromules
	      case 0: NBRANCH = 0; break;
		// single stromules
	      case 1: NBRANCH = 2;
		theta = (dist == 0 && sample == 0 ? 0 : RND*3.14159*2);
		x1[0] = x10; y1[0] = y10; x2[0] = x10+(r+L)*cos(theta); y2[0] = y10+(r+L)*sin(theta);
		theta = RND*3.14159*2;
		x1[1] = x20; y1[1] = y20; x2[1] = x20+(r+L)*cos(theta); y2[1] = y20+(r+L)*sin(theta);
		break;
		// branched stromules
	      case 2: NBRANCH = 6;
		theta = (dist == 0 && sample == 0 ? 0 : RND*3.14159*2);
		x1[0] = x10; y1[0] = y10; x2[0] = x10+(r+L/3)*cos(theta); y2[0] = y10+(r+L/3)*sin(theta);
		x1[1] = x2[0]; y1[1] = y2[0]; x2[1] = x1[1]+(L/3)*cos(theta+PI+2*PI/3); y2[1] = y1[1]+(L/3)*sin(theta+PI+2*PI/3);
		x1[2] = x2[0]; y1[2] = y2[0]; x2[2] = x1[2]+(L/3)*cos(theta+PI-2*PI/3); y2[2] = y1[2]+(L/3)*sin(theta+PI-2*PI/3);
		theta = RND*3.14159*2;
		x1[3] = x20; y1[3] = y20; x2[3] = x20+(r+L/3)*cos(theta); y2[3] = y20+(r+L/3)*sin(theta);
		x1[4] = x2[3]; y1[4] = y2[3]; x2[4] = x1[4]+(L/3)*cos(theta+PI+2*PI/3); y2[4] = y1[4]+(L/3)*sin(theta+PI+2*PI/3);
		x1[5] = x2[3]; y1[5] = y2[3]; x2[5] = x1[5]+(L/3)*cos(theta+PI-2*PI/3); y2[5] = y1[5]+(L/3)*sin(theta+PI-2*PI/3);
		break;
		// double stromules
	      case 3: NBRANCH = 4;
		theta = (dist == 0 && sample == 0 ? 0 : RND*3.14159*2);
		x1[0] = x10; y1[0] = y10; x2[0] = x10+(r+L/2)*cos(theta); y2[0] = y10+(r+L/2)*sin(theta);
		x1[1] = x10; y1[1] = y10; x2[1] = x10+(r+L/2)*cos(theta+PI); y2[1] = y10+(r+L/2)*sin(theta+PI);
		theta = RND*3.14159*2;
		x1[2] = x20; y1[2] = y20; x2[2] = x20+(r+L/2)*cos(theta); y2[2] = y20+(r+L/2)*sin(theta);
		x1[3] = x20; y1[3] = y20; x2[3] = x20+(r+L/2)*cos(theta+PI); y2[3] = y20+(r+L/2)*sin(theta+PI);
		break;
	      }

	      if(dist == 0) NBRANCH /= 2;
	      
     	      if(sample != -1)
		{
		  fp = fopen("geom-snaps.csv", "a");
		  for(i = 0; i < NBRANCH; i++)
		    {
		      fprintf(fp, "%i,%f,%i,%f,%f,%f,%f\n", expt, dist, sample, x1[i], y1[i], x2[i], y2[i]);
		    }
		  fclose(fp);
		}

	      if(sample == 0)
		{
		  sprintf(fstr, "cellmodelmany-%i-%.0f.txt", expt, dist);
	          fp = fopen(fstr, "w");
		}
	      for(cx = 0; cx < DIM; cx++)
		{
		  for(cy = 0; cy < DIM; cy++)
		    {
		      x = cx-DIM/2; y = cy-DIM/2;

		      if((x-x10)*(x-x10)+(y-y10)*(y-y10) <= r*r || (x-x20)*(x-x20)+(y-y20)*(y-y20) <= r*r) IR[cy*DIM+cx] = PA[cy*DIM+cx] = 0;
		      else
			{
			  // distance to first main plastid XXX FIX IN ORIGINAL
			  d = sqrt((x-x10)*(x-x10)+(y-y10)*(y-y10)) - r;
			  // distance to second main plastid
			  thisd = sqrt((x-x20)*(x-x20)+(y-y20)*(y-y20)) - r;
			  if(thisd < d) d = thisd;
			  // check distances to any stromule segments to see if these are shorter
			  for(i = 0; i < NBRANCH; i++)
			    {
			      thisd = pDistance(x, y, x1[i], y1[i], x2[i], y2[i]);
			      if(thisd < d) d = thisd;
			    }
			  IR[cy*DIM+cx] = d;

			  countcollision = 0;
			  for(theta = 0; theta < 2*3.1415; theta += dtheta)
			    {
			      thiscollision = 0;
			      for(i = 0; i < NBRANCH && thiscollision == 0; i++)
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
				  // small segment for plastid 1
				  tx1 = x10 + r*sin(theta2); ty1 = y10 + r*cos(theta2);
				  tx2 = tx1 + dtheta2*r*sin(theta2-PI/2); ty2 = ty1 + dtheta2*r*cos(theta2-PI/2);
				  // detect collision
				  collision = get_line_intersection(x, y, x+2*DIM*sin(theta), y+2*DIM*cos(theta), tx1, ty1, tx2, ty2, &ix, &iy);
				  if(collision)
				    {
				      // if collision point lies within RHO, record it
				      if((x-ix)*(x-ix) + (y-iy)*(y-iy) < RHO*RHO)
					{
					  thiscollision = 1;
					  break;
					}
				    }
				  if(dist != 0)
				    {
				      // small segment for plastid 2
				      tx1 = x20 + r*sin(theta2); ty1 = y20 + r*cos(theta2);
				      tx2 = tx1 + dtheta2*r*sin(theta2-PI/2); ty2 = ty1 + dtheta2*r*cos(theta2-PI/2);
				      // detect collision
				      collision = get_line_intersection(x, y, x+2*DIM*sin(theta), y+2*DIM*cos(theta), tx1, ty1, tx2, ty2, &ix, &iy);
				      if(collision)
					{
					  // if collision point lies within RHO, record it
					  if((x-ix)*(x-ix) + (y-iy)*(y-iy) < RHO*RHO)
					    {
					      thiscollision = 1;
					      break;
					    }
					}
				    }

				}
			      countcollision += thiscollision;
			    }
			  PA[cy*DIM+cx] = (double)countcollision/(theta/dtheta);
			}
		      histIR[(int)(IR[cy*DIM+cx]/(2*DIM)*BINS)]++;
		      histPA[(int)(PA[cy*DIM+cx]*BINS)]++;
		      if(sample == 0)
			{
			  fprintf(fp, "%i %i %f %f\n", cx, cy, IR[cy*DIM+cx], PA[cy*DIM+cx]);
			}
		    }
		}
	      if(sample == 0)
		{
		  fclose(fp);
		}
	    }
	
	  fp = fopen("hists-many.csv", "a");
	  for(i = 0; i <= BINS; i++)
	    fprintf(fp, "%i,%f,%i,%f,%f\n", expt, dist, 1, (double)i/BINS, histIR[i]/nsample);
	  for(i = 0; i <= BINS; i++)
	    fprintf(fp, "%i,%f,%i,%f,%f\n", expt, dist, 2, (double)i/BINS, histPA[i]/nsample);
	  fclose(fp);
	}
    }
  
  return 0;
}
  
		  
