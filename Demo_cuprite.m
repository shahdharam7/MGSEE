%% demo_MGSEE
clc;
close all;
clear all;

%% Image Read
s=load('Cuprite.mat');  % link for data source : https://rslab.ut.ac.ir/data
p=s.nRow;
q=s.nCol;
Bands=188;
Y=s.Y;
x=hyperConvert3d(Y,p,q,Bands);

%% Virtual Dimension
VD=12;

%% MGSEE algorithm
[endmemberindex] = MGSEE(x,VD);
endmemberindex_CPM=change_index(endmemberindex,p,q);

%% VCA algorithm
[U_VCA,e_index,snrEstimate]=hyperVca(Y,VD);
endmemberindex_VCA=change_index(e_index,p,q);    

%% gt compare
t1=load('groundTruth_Cuprite_nEnd12.mat');
gt=t1.M;
n1=gt(3:103,:);
n2=gt(114:147,:);
n3=gt(168:220,:);
gt=[n1;n2;n3];
[gt_m,gt_n]=size(gt);

for i=1:gt_n
    for j=1:Bands
        extracted_VCA(j,i)=x(endmemberindex_VCA(i,1),endmemberindex_VCA(i,2),j);
        extracted_CPM(j,i)=x(endmemberindex_CPM(i,1),endmemberindex_CPM(i,2),j);
    end
end

%% SAM Calculation
ex_n=gt_n;
store_VCA=[0,0];
store_CPM=[0,0];
sam_VCA=0;
sam_CPM=0;
sam_total_VCA=0;
sam_total_CPM=0;

for i=1:gt_n
    for j=1:ex_n
        Mat_SAM_CPM(i,j)=real(acos(dot(gt(:,i),extracted_CPM(:,j))/(norm(gt(:,i)*norm(extracted_CPM(:,j))))));
        Mat_SAM_VCA(i,j)=real(acos(dot(gt(:,i),extracted_VCA(:,j))/(norm(gt(:,i)*norm(extracted_VCA(:,j))))));
    end
end

for i=1:gt_n
    %CPM
    [max_value1,mrow]=min(Mat_SAM_CPM);
    [max_value,col_CPM]=min(max_value1);
    sam_total_CPM=sam_total_CPM+max_value;
    sam_CPM=[sam_CPM;max_value];
    row_CPM=mrow(col_CPM);
    s1=[row_CPM,col_CPM];
    store_CPM=[store_CPM;s1];
    save_CPM(row_CPM)=max_value;
    Mat_SAM_CPM(row_CPM,:)=[100*ones];
    Mat_SAM_CPM(:,col_CPM)=[100*ones];
    %VCA
    [max_value1,mrow]=min(Mat_SAM_VCA);
    [max_value,col_VCA]=min(max_value1);
    sam_total_VCA=sam_total_VCA+max_value;
    sam_VCA=[sam_VCA;max_value];
    row_VCA=mrow(col_VCA);
    s1=[row_VCA,col_VCA];
    store_VCA=[store_VCA;s1];
    save_VCA(row_VCA)=max_value;
    Mat_SAM_VCA(row_VCA,:)=[100*ones];
    Mat_SAM_VCA(:,col_VCA)=[100*ones];
end

rms_sae=[rms(save_CPM);
    rms(save_VCA)];
rms_sae = radtodeg(rms_sae);

disp('RMSSAE of VCA');
disp(rms_sae(2));
disp('RMSSAE of CPM');
disp(rms_sae(1));