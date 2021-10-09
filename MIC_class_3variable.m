function [MIC,mycertain,certain_seg,bestc]=MIC_class_3variable(data,B,c,n)
% n is the numbers of sample
% 1th column of data is y,2th and 3th column are x1 and x2 respectively
% the entries of y only contain 1 and 2
% mutual_I_2 represent fixed the first column of data(x1)
% mutual_I_1 represent fixed the second column of data(x2)
% n^B=n^0.6
% c is max parameter of clumps (c*x)
% "mycertain" represent which x is equipartition,1 is the first x,2 is the
% second x; "certain_seg" represent how many segments are divided of the
% certained x; "bestc" the segments of another x.
%
class_num=length(unique(data(:,1)));
% multi=1;
multi=class_num;
% multi=log2(class_num);
max_seg=round((n/2)^B/multi);
class_num=int32(class_num);
if max_seg<2
    max_seg=2;
end
mutual_I_2=zeros(max_seg-1,max_seg-1);
mutual_I_1=zeros(max_seg-1,max_seg-1);
last_c2=zeros(max_seg+1,max_seg+1);
last_c1=zeros(max_seg+1,max_seg+1);
sample_n=n;
n=int32(n);
% randnum=randperm(n)';
% data=data(randnum,:);
% [~,pos1]=sortrows(data,2);
% [~,pos2]=sortrows(data,3);
[value1,pos1]=sortrows(data,2);
[value2,pos2]=sortrows(data,3);
D1=nan(n,2);
D1=int32(D1);
vector_y_x1=data(pos2,1);
vector_y_x1=int32(vector_y_x1);
vector_y_x2=data(pos1,1);
vector_y_x2=int32(vector_y_x2);
for i=2:max_seg
    %     Q_x=equipartitionYaxis(int32(i),n);
    Q_x1=equipartitionYaxis2c(value1(:,2),int32(i),n);
    Q_x2=equipartitionYaxis2c(value2(:,3),int32(i),n);
    avg=int32(sample_n/(round((sample_n/i)^B/multi)*c));
    sub_max_seg=int32(round((sample_n/i)^B/multi));
    if max(Q_x1)==i
        D1(pos1,1)=Q_x1;
        vector1=D1(pos2,1);
        best_c1=getsuper3var(vector1,vector_y_x1,avg,n,class_num);
        len1=int32(length(best_c1)-2);
        [temI2,tem_c2]=getmutualI3var(vector1,vector_y_x1,best_c1,sub_max_seg,int32(i),n,len1,class_num);
        mutual_I_2(i-1,1:length(temI2))=temI2';
        last_c2(i-1,1:length(tem_c2))=tem_c2;
    end
    %     avg=int32(n/(round(B/i)*c));
    %     best_c1=getsuper(vector1,vector_y_x1,avg,n);
    %     best_c2=getsuper(vector2,vector_y_x2,avg,n);
    %     [temI2,tem_c2]=getmutualI(vector1,vector_y_x1,best_c1,sub_max_seg,int32(i),n,len1);
    %     [temI1,tem_c1]=getmutualI(vector2,vector_y_x2,best_c2,sub_max_seg,int32(i),n,len2);
    if max(Q_x2)==i
        D1(pos2,2)=Q_x2;
        vector2=D1(pos1,2);
        best_c2=getsuper3var(vector2,vector_y_x2,avg,n,class_num);
        len2=int32(length(best_c2)-2);
        [temI1,tem_c1]=getmutualI3var(vector2,vector_y_x2,best_c2,sub_max_seg,int32(i),n,len2,class_num);
        mutual_I_1(1:length(temI1),i-1)=temI1;
        last_c1(1:length(tem_c1),i-1)=tem_c1';
    end
end
MIC_min=min(min([mutual_I_1,mutual_I_2]));
MIC_max=max(max([mutual_I_1,mutual_I_2]));
if abs(MIC_max)>abs(MIC_min)
    MIC=MIC_max;
else
    MIC=MIC_min;
end
position=find(mutual_I_2==MIC, 1);
if ~isempty(position)
    pos_row=mod(position(1)-1,round((sample_n/2)^B/multi)-1)+1;
    pos_col=floor((position(1)-1)/(round((sample_n/2)^B/multi)-1))+1;
    mycertain=1;
    certain_seg=pos_row+1;
    bestc=last_c2(pos_row,2:pos_col+1);
else
    position=find(mutual_I_1==MIC, 1);
    pos_row=mod(position(1)-1,round((sample_n/2)^B/multi)-1)+1;
    pos_col=floor((position(1)-1)/(round((sample_n/2)^B/multi)-1))+1;
    mycertain=2;
    certain_seg=pos_col+1;
    bestc=last_c1(2:pos_row+1,pos_col)';
end
end
