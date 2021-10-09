#include "mex.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#define NULL 0
#define LEN sizeof(struct clumb)
struct clumb
{int subc;
 struct clumb *next;
};
struct clumb *insert(struct clumb *head,struct clumb *stud)
{
 struct clumb *sub_node0,*sub_node1,*sub_node2;  
 sub_node1=head;
 sub_node0=stud;
 sub_node2=NULL;
 if(head==NULL){
     head=sub_node0;
     sub_node0->next=NULL;
 }else{
     while((sub_node0->subc>sub_node1->subc)&&(sub_node1->next!=NULL)){
         sub_node2=sub_node1;
         sub_node1=sub_node1->next;         
     }
     if(sub_node0->subc<=sub_node1->subc){
         if(head==sub_node1) head=sub_node0;
         else sub_node2->next=sub_node0;
         sub_node0->next=sub_node1;   
     }else{
         sub_node1->next=sub_node0;
         sub_node0->next=NULL;    
     }
 }
 return(head);
}
struct clumb *del(struct clumb *head,int num)
{
 struct clumb *sub_node1,*sub_node2;
 sub_node1=head;
//  sub_node1=sub_node1->next;
 if(head==NULL) mexPrintf("list is null\n");
 else
 {while(num!=sub_node1->subc && sub_node1->next!=NULL){
     sub_node2=sub_node1;
     sub_node1=sub_node1->next;  
 }
 if(num==sub_node1->subc){
     if(sub_node1==head) head=sub_node1->next;
     else sub_node2->next=sub_node1->next;
     mxFree(sub_node1);
 }
 else mexPrintf("not been found\n");
}
 return(head);
}
void print(struct clumb *head)
{
 struct clumb *node;
 node=head;
 if(head!=NULL)
     do
     {
     mexPrintf("%d\n",node->subc);
     node=node->next;
     } while(node!=NULL);
}
void release(struct clumb *head)
{
 struct clumb *node;
//  node=head;
 if(head!=NULL)
     do
     {
     node=head;
     head=head->next;
     mxFree(node);
     } while(head!=NULL);
}
void listtomatric(struct clumb *head,int *submat,int len_mat)
{
 struct clumb *node;
 int subcount=0;
 node=head;
 if(head!=NULL)
     do
     {
     if(subcount!=0) submat[subcount-1]=node->subc;
     subcount=subcount+1;
     node=node->next;
     } while(node!=NULL);
}
double myentropy(int mytable[],int tab_len,int num_sample)
{
    int i;
    double entr=0.0;
    for(i=0;i<tab_len;i++){
        if(mytable[i]!=0){
            entr=entr+(((double)mytable[i])/((double)num_sample))*(log(((double)mytable[i])/((double)num_sample))/log(2)); 
        }
    }
    entr=entr*(-1);
    return entr;
}
double mutual_I(int vector_x[],struct clumb *head,int certain,int num_sample,int vector_y[],int len_seg,int class_num)
{
    int i,j,jk,subk=0,*subcertain,*class_,*temseg,*array_cer_y,*array_cer,*array_seg_y,*array_seg;
    int *array_seg_cer,*array_seg_cer_y,*array_y;
    double myI,H1,H2,H3,H4,H5,H6,H7,maxH;
    array_cer_y=(int *) mxCalloc(class_num*certain,sizeof(int));
    array_cer=(int *) mxCalloc(certain,sizeof(int));
    subcertain=(int *) mxCalloc(certain,sizeof(int));
    class_=(int *) mxCalloc(class_num,sizeof(int));
    array_seg_y=(int *) mxCalloc(class_num*len_seg,sizeof(int));
    array_seg=(int *) mxCalloc(len_seg,sizeof(int));
    array_seg_cer=(int *) mxCalloc(certain*len_seg,sizeof(int));
    array_seg_cer_y=(int *) mxCalloc(class_num*certain*len_seg,sizeof(int));
    array_y=(int *) mxCalloc(class_num,sizeof(int));
    for(i=0;i<class_num*certain*len_seg;i++){
        if(i<class_num) array_y[i]=0;
        if(i<certain){
            array_cer[i]=0;
            subcertain[i]=i+1;
        }
        if(i<class_num) class_[i]=i+1;
        if(i<class_num*certain) array_cer_y[i]=0;
        if(i<class_num*len_seg) array_seg_y[i]=0;
        if(i<len_seg) array_seg[i]=0;
        if(i<certain*len_seg) array_seg_cer[i]=0;
        array_seg_cer_y[i]=0;
    }
    temseg=(int *) mxCalloc(len_seg,sizeof(int));
    listtomatric(head,temseg,len_seg);
    for(i=0;i<num_sample;i++){
        for(j=0;j<certain;j++){
            if(vector_x[i]==subcertain[j]){
                array_cer[j]=array_cer[j]+1;
                if(i<temseg[subk]){
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                    array_seg[subk]=array_seg[subk]+1;
                }else{
                    subk=subk+1;
                    array_seg[subk]=array_seg[subk]+1;
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                }
                for (jk=0;jk<class_num;jk++){
                    if(vector_y[i]==class_[jk]){
                        array_seg_y[jk+subk*class_num]=array_seg_y[jk+subk*class_num]+1;
                        array_y[jk]=array_y[jk]+1;
                        array_cer_y[jk+class_num*j]=array_cer_y[jk+class_num*j]+1;
                        array_seg_cer_y[jk+class_num*j+subk*class_num*certain]=array_seg_cer_y[jk+class_num*j+subk*class_num*certain]+1; 
                        break;
                    }
                }
              break;
            }
        }
    }
    H1=myentropy(array_seg_cer_y,class_num*certain*len_seg,num_sample); 
    H2=myentropy(array_seg_y,class_num*len_seg,num_sample);
    H3=myentropy(array_cer_y,class_num*certain,num_sample);
    H4=myentropy(array_seg_cer,certain*len_seg,num_sample);
    H5=myentropy(array_seg,len_seg,num_sample);
    H6=myentropy(array_cer,certain,num_sample);
    H7=myentropy(array_y,class_num,num_sample);
//     maxH=len_seg;
//     if (maxH<certain){
//         maxH=certain;
//     }
//     if (maxH<class_num){
//         maxH=class_num;
//     }
//     myI=2*(H2+H3+H4-H1-H5-H6-H7)/(log(len_seg)/log(2)+log(certain)/log(2)+log(class_num)/log(2)-log(maxH)/log(2));
//     maxH=len_seg;
//     if (maxH>class_num){
//         maxH=class_num;
//     }
    maxH=class_num;
    myI=(H2+H3+H4-H1-H5-H6-H7)/(log(maxH)/log(2));
    return myI;
    mxFree(array_cer_y);
    mxFree(array_cer);
    mxFree(array_seg_y);
    mxFree(array_seg);
    mxFree(class_);
    mxFree(array_seg_cer_y);
    mxFree(array_seg_cer);
    mxFree(array_y);
    mxFree(subcertain);
    mxFree(temseg);
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int *vector_x, *vector_y,*best_c,*out2,*outtem;
    int num_sample,segm,certain,len_bestc,class_num,len_segm,i,j,count,num;
    double *myI,I0,tem_I;
    struct clumb *head_mybestc,*head_tembestc,*node1,*node2,*node_tem;
    vector_x=mxGetData(prhs[0]);
    vector_y=mxGetData(prhs[1]);
    best_c=mxGetData(prhs[2]);
    segm=*((int *)mxGetData(prhs[3]));
    certain=*((int *)mxGetData(prhs[4]));
    num_sample=*((int *)mxGetData(prhs[5]));
    len_bestc=*((int *)mxGetData(prhs[6]));
    class_num=*((int *)mxGetData(prhs[7]));
    if(len_bestc>segm-1) len_segm=segm-1;
    else len_segm=len_bestc;
    plhs[0]=mxCreateDoubleMatrix(len_segm,1,mxREAL);
    myI=mxGetPr(plhs[0]);
    plhs[1]=mxCreateNumericMatrix(len_segm+2,1,mxINT32_CLASS,mxREAL);
    out2=mxGetData(plhs[1]);
    head_mybestc=NULL;
    head_tembestc=NULL;
    node1=node2=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[0];
    head_mybestc=node1;
    node2=node1;
    node1=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[len_bestc+1];
    node2->next=node1;
    node1->next=NULL;
//     
    node1=node2=(struct clumb *) mxMalloc(LEN);
    for(i=1;i<len_bestc+1;i++){
        node1->subc=best_c[i];
        if(i==1) head_tembestc=node1;
        else node2->next=node1;
        node2=node1;
        node1=(struct clumb *) mxMalloc(LEN);
    }
    node2->next=NULL;     
    count=0;
    out2[0]=0;
    for(i=0;i<len_segm;i++){
        I0=0.0;
        node_tem=head_tembestc;
        if(head_tembestc!=NULL)
        do
        {
        j=node_tem->subc;
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=j;
        head_mybestc=insert(head_mybestc,node1);
        tem_I=mutual_I(vector_x,head_mybestc,certain,num_sample,vector_y,count+2,class_num);
        if (fabs(tem_I)>fabs(I0)){
            I0=tem_I;
            num=j;
        }
        head_mybestc=del(head_mybestc,j);
        node_tem=node_tem->next;
        }while(node_tem!=NULL);
        myI[count]=I0;
        count++;    
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=num;
        out2[i+1]=num;
        head_mybestc=insert(head_mybestc,node1);
        head_tembestc=del(head_tembestc,num);
    }
    out2[len_segm+1]=num_sample;
    release(head_mybestc);
    release(head_tembestc);
    return;
}