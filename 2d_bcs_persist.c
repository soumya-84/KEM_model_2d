#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ran2.c"
#define L 24
float expd(void);
long sd=-972434472;
main()
{
  int i,j,k,count2,iw,count,n,NRUN=10,t,TMAX=10000,ii,jj,ni,nj,id,xl,yl,**group,**super_spin,M=7,tt,flag,ratio,RELAX=10000,counti,countj,delta,idelta,inorm,LARGE=L*L,clg_list[51],nr,a,used[51],limit,before,after,**stat,**flip,count1,count0,flip_number,en,absorb=0,en_h,flag1;
  float sum,r,p,mu,*op,tot_sum,*s_op,*time_avg,zero,tot_flip,tot_zero,*norm,*prob_del,factor[M][M],total_factor=0.,tot_op,new_sum,*per,*per1,*per0,*domain,*activity,*super_active,*abs_dist;
  FILE *fpt1,*fpt2,*fpt3;

 ////////////////////////////can change NRUN depending on run-time


  ratio=L/M;
 //   norm=(float *)calloc(LARGE, sizeof(float));
//  prob_del=(float *)calloc(LARGE, sizeof(float));
  
//  fpt3=fopen("elec_clg_list_3","r");
//  fpt1=fopen("mag_L30_p0.12n","w");
//  fpt2=fopen("flip_prob_L30_p0.12n","w");
 
///////   p=0.5;

     stat = (int **)calloc(L, sizeof(int*));
     group = (int **)calloc(L,sizeof(int*));
     flip = (int **)calloc(L,sizeof(int*));
     op = (float *)calloc(TMAX,sizeof(float));
     s_op=(float *)calloc(TMAX,sizeof(float));
     time_avg=(float *)calloc(TMAX,sizeof(float));
     per=(float *)calloc(TMAX,sizeof(float));
     per1=(float *)calloc(TMAX,sizeof(float));
     per0=(float *)calloc(TMAX,sizeof(float));
     domain=(float *)calloc(TMAX,sizeof(float));
     activity=(float *)calloc(TMAX,sizeof(float));
     super_active =(float *)calloc(TMAX,sizeof(float));
     abs_dist = (float *)calloc(TMAX,sizeof(float));
 //    super_spin = (int **)calloc(M,sizeof(int*));


//  for(i=0;i<51;i++) {fscanf(fpt3,"%d %d\n",&clg_list[i],&a); }

  
//  limit=50;



  for(i=0;i<L;i++) 
  {
     stat[i]=(int *)calloc(L, sizeof(int));
     group[i]=(int *)calloc(L,sizeof(int));
     flip[i]=(int *)calloc(L,sizeof(int));
  }

// for(p=0.;p<.2;p=p+0.02)             ///////////////////////////////  the range of p can be changed here/////////////////////////
// for(limit=1;limit<40;limit=limit+1)
 {
  tot_op=0;

  p=0.;
  

  for(i=0;i<TMAX;i++) {op[i]=0; s_op[i]=0; time_avg[i]=0; per[i]=0; per1[i]=0; per0[i]=0; domain[i]=0; activity[i]=0; super_active[i]=0; abs_dist[i]=0;}
 
/////////////relaxation//////////////////////////
/*
     for(i=0;i<L;i++)
     {
       for(j=0;j<L;j++)
       {
          if(ran2(&sd)<.5) {stat[i][j]=1;}
          else {stat[i][j]=-1;}
       }
     } 

    for(t=0;t<RELAX;t++)
    {
        for(j=0;j<L*L;j++)
         {
            ii=ran2(&sd)*L; jj=ran2(&sd)*L;
            r=ran2(&sd);
            if(r<0.25) {ni=(ii-1+L)%L; nj=jj;}
            else if(r<0.5) {ni=(ii+1)%L; nj=jj;}
            else if(r<0.75) {ni=ii; nj=(jj-1+L)%L;}
            else {ni=ii; nj=(jj+1)%L;}


         //   ni=ran2(&sd)*L; nj=ran2(&sd)*L;    //mean field limit

            if(ran2(&sd)<p) {mu=-1;} else {mu=+1;}

           if(matrix[ii][jj]*matrix[ni][nj]==1) 
           { 
              stat[ii][jj]=stat[ii][jj]+mu*stat[ni][nj];
              if(stat[ii][jj]<-1) {stat[ii][jj]=-1;}
              if(stat[ii][jj]>1) {stat[ii][jj]=1;}
           }
         }
    }
*/
/////////////////////////////////////////////////

  zero=0; tot_flip=0; flag=0; tot_zero=0;

  for(n=0;n<NRUN;n++)
  {
     for(i=0;i<L;i++)
     {
       for(j=0;j<L;j++)
       {
          r=ran2(&sd);
          if(r<.33) {stat[i][j]=1;}
          else if(r<.66){stat[i][j]=-1;}
          else {stat[i][j]=0;}
          group[i][j]=0;
          flip[i][j]=1;
       }
     }   
     tt=0;
     flag1=1;
     for(t=0;t<TMAX;t++)
     {
         flip_number=0;
         for(j=0;j<L*L;j++)
         {
            ii=ran2(&sd)*L; jj=ran2(&sd)*L;

            r=ran2(&sd);
            if(r<0.25) {ni=(ii-1+L)%L; nj=jj;}
            else if(r<.5) {ni=(ii+1)%L; nj=jj;}
            else if(r<0.75) {ni=ii; nj=(jj-1+L)%L;}
            else {ni=ii; nj=(jj+1)%L;}


       //      ni=ran2(&sd)*L; nj=ran2(&sd)*L;    //mean field limit

          if(ran2(&sd)<p) {mu=-1;} else {mu=1;}
    //      if(group[ii][jj]<=limit)
          { 
            before=stat[ii][jj];
            stat[ii][jj]=stat[ii][jj]+mu*stat[ni][nj];
            if(stat[ii][jj]<-1) {stat[ii][jj]=-1;}
            if(stat[ii][jj]>1) {stat[ii][jj]=1;}
            after=stat[ii][jj];
            if(before!=after) {group[ii][jj]++; flip[ii][jj]=0; flip_number++;}
          }
         }

         sum=0;

        ii=0; jj=0;
         zero=0; counti=0; count0=0; count1=0; en=0; en_h=0;
         for(i=0;i<L;i++)
         {
       
            for(j=0;j<L;j++)
            {
                
              //   printf("i=%d j=%d ii=%d jj=%d\n",i,j,ii,jj);
                sum+=stat[i][j];
                counti+=flip[i][j];
                if(flip[i][j]==1 && stat[i][j]==1) {count1++;}
                if(flip[i][j]==1 && stat[i][j]==0) {count0++;}


                if(stat[i][j]!=stat[(i+1)%L][j]) {en++;}
                if(stat[i][j]!=stat[i][(j+1)%L]) {en++;}

                if(stat[i][j]*stat[(i+1)%L][j]==-1) {en_h++;}
                if(stat[i][j]*stat[i][(j+1)%L]==-1) {en_h++;}


            }
          
         }
         
         op[t]+=fabs(sum)/(L*L);
         per[t]+=(float)counti/(L*L);
         per1[t]+=(float)count1/(L*L);
         per0[t]+=(float)count0/(L*L);
         domain[t]+=(float)en/(L*L);
         activity[t]+=(float)flip_number/(L*L);
         super_active[t]+=(float)en_h/(L*L);

         if(flag1==1 && fabs(sum)==L*L)
         {
            flag1=0;
            abs_dist[t]++; //printf("t=%d n=%d\n",t,n);
         /////   t=TMAX;      //coming out of time loop
         }


     }
//end of time loop

//  tot_flip+=(float)flag;
  
    if((float)sum/(L*L)<0.001) {absorb++;}



  }
//end of realisation

 //// printf("%d %f\n",limit,op[TMAX-1]/NRUN);

 //   printf("%f %f %f %f\n",p,tot_zero/(TMAX*NRUN),(float)flag/(NRUN*TMAX),tot_op/(NRUN*TMAX));
 }
//end of p loop
   
  for(i=0;i<TMAX;i++)
  {
   ///////// printf("%d %f\n",i,op[i]);
///   if(abs_dist[i]>0) { printf("%d %f\n",i,abs_dist[i]/NRUN);}
     printf("%d %f %f %f %f %f %f %f %f %f\n",i,op[i]/NRUN,per[i]/NRUN,per1[i]/NRUN,per0[i]/NRUN,domain[i]/NRUN,activity[i]/NRUN,super_active[i]/NRUN,(float)absorb/NRUN,abs_dist[i]/NRUN);
  //   fprintf(fpt2,"%d %f\n",i,time_avg[i]/NRUN);
  }



  for(i=0;i<L;i++)
  {
    for(j=0;j<L;j++)
    {
   ///    printf("%d ",stat[i][j]);
    }
  ///  printf("\n");
  }

}


float expd(void)
{

   float lam, x,y,z;
  
   lam=0.01; 

   z=ran2(&sd);
   x=log(1-z)/(-1*lam);
 
   return(x);
}







