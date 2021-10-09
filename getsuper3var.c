#include "mex.h"
#include "math.h"
void find_index(int y[],int mx,int rows,int len[],int *ponit_pos)
{
    int i,j,n,m,count;
    for (j=0;j<rows;j++){
        n=0;
        m=len[j]-1;
        i=(n+m)/2;
        count=1;
        while(n!=m && y[j+rows*i]!=mx && count<len[j]){
            if(y[j+rows*i]<mx) n=i+1;
            if(y[j+rows*i]>mx) m=i-1;
            i=(n+m)/2;
            count++;
        }
        if(y[j+rows*i]==mx){
            ponit_pos[0]=i+1;
            ponit_pos[1]=j;
            break;
        }
    }
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mxArray *vector_x_class, *rank_x;
    int *vector_x, *vector_y, *ptr_vector_x_class, *ptr_rank_x,*best_c;
    int *out1,*rank_,rank_sample,rank_c,num_sample,avg;
    int num1,num2,temc1,*out4,class_num;
    int *point_pos,i,j,k=1,class_sample,*len_sample,*len_sample2,max_num;
    int if_fix=1;
    rank_sample=1;
    vector_x=mxGetData(prhs[0]);
    vector_y=mxGetData(prhs[1]);
    avg=*((int *)mxGetData(prhs[2]));
    num_sample=*((int *)mxGetData(prhs[3]));
    class_num=*((int *)mxGetData(prhs[4]));
    rank_=(int *) mxCalloc(num_sample,sizeof(int));
    best_c=(int *) mxCalloc(num_sample,sizeof(int));
    len_sample=(int *) mxCalloc(class_num,sizeof(int));
    len_sample2=(int *) mxCalloc(class_num,sizeof(int));
    point_pos=(int *) mxCalloc(2,sizeof(int));
    for (i=0;i < class_num; i++){
        len_sample[i]=0;
        len_sample2[i]=0;
    }
    for (i = 0; i < num_sample; i++){
        rank_[i]=i+1;
        best_c[i]=0;       
        for (j=0; j < class_num; j++){
            if(vector_y[i]==j+1){
                len_sample[j]++;
                break;
            }
        }
    }
    max_num=0;
    for (i=0; i < class_num; i++){
        if (max_num<len_sample[i]){
            max_num=len_sample[i];
        }
    }
    vector_x_class=mxCreateNumericMatrix(class_num,max_num,mxINT32_CLASS,mxREAL);
    rank_x=mxCreateNumericMatrix(class_num,max_num,mxINT32_CLASS,mxREAL);
    ptr_vector_x_class=mxGetData(vector_x_class);
    ptr_rank_x=mxGetData(rank_x);
    for (i = 0; i < num_sample; i++){
        for (j = 0; j < class_num; j++){
            if(vector_y[i]==j+1){
                ptr_vector_x_class[j+len_sample2[j]*class_num]=vector_x[i];
                ptr_rank_x[j+len_sample2[j]*class_num]=rank_[i];
                len_sample2[j]++;
                break;
            }
        }
    }
    point_pos[0]=0;
    point_pos[1]=0;
    find_index(ptr_rank_x,rank_sample,class_num,len_sample,point_pos);
    rank_c=ptr_rank_x[point_pos[1]+(point_pos[0]-1)*class_num];
    while(rank_sample<num_sample){
        if ((point_pos[0] < len_sample[point_pos[1]]) && ptr_vector_x_class[point_pos[1]+(point_pos[0]-1)*class_num]==ptr_vector_x_class[point_pos[1]+point_pos[0]*class_num]){        
            point_pos[0]++;
            rank_c=ptr_rank_x[point_pos[1]+(point_pos[0]-1)*class_num];
            rank_sample=rank_c+1;
        }else{
            rank_sample=rank_c+1;
            find_index(ptr_rank_x,rank_sample,class_num,len_sample,point_pos);
            if (if_fix==1){
                if ((rank_c-best_c[k-1]) < avg){
                   num1= rank_c-best_c[k-1];
                   temc1=rank_c;
                   if_fix=0;  
                }else{
                    best_c[k]=rank_c;
                    k=k+1;
                    if_fix=1;
                }
            }else{
                num2=rank_c-temc1;
                if ((num2+num1) < avg){
                    num1=num2+num1;
                    temc1=rank_c;
                    if_fix=0;
                }else{ 
                    if ((avg-num1)>(num1+num2)-avg){
                        best_c[k]=rank_c;
                        if_fix=1;
                        k=k+1;
                    }else{
                        best_c[k]=temc1;
                        k=k+1;
                        if (num2>=avg){
                            best_c[k]=rank_c;
                            k=k+1;
                            if_fix=1;
                        }else{
                            num1=num2;
                            temc1=rank_c;
                            if_fix=0;
                        }
                    }
                }
            }
            rank_c=ptr_rank_x[point_pos[1]+(point_pos[0]-1)*class_num];  
         }
    }
    if (best_c[k-1]!=num_sample){
        best_c[k]=num_sample;
        plhs[0]=mxCreateNumericMatrix(1,k+1,mxINT32_CLASS,mxREAL);
        out1=mxGetData(plhs[0]);
        for (i=0;i<k+1;i++){
              out1[i]=best_c[i];
        }
    }else{
        plhs[0]=mxCreateNumericMatrix(1,k,mxINT32_CLASS,mxREAL);
        out1=mxGetData(plhs[0]);
        for (i=0;i<k;i++){
              out1[i]=best_c[i];
        }
    }
    mxFree(rank_);
    mxFree(best_c);
    mxFree(len_sample);
    mxFree(len_sample2);
    mxFree(point_pos);
    mxDestroyArray(vector_x_class);
    mxDestroyArray(rank_x);
    return;
}