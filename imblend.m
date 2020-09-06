function  output = imblend( source, mask, target, transparent )
%Source, mask, and target are the same size (as long as you do not remove
%the call to fiximages.m). You may want to use a flag for whether or not to
%treat the source object as 'transparent' (e.g. taking the max gradient
%rather than the source gradient).

%output = source .* mask + target .* ~mask;

T_size=size(target);
M_size=size(mask);
S_size=size(source);
A_size=T_size(1)*T_size(2);


% Solve the problem of not having neighbours in the edges
Maskedge=0;
if( (any(mask(1,:,1)==1)) || (any(mask(:,1,1)==1)) || (any(mask(M_size(1),:,1)==1)) || (any(mask(:,M_size(2),1)==1)) )
    
    temp_target=zeros(T_size(1)+2,T_size(2)+2,3);
    temp_mask=zeros(M_size(1)+2,M_size(2)+2,3);
    temp_source=zeros(S_size(1)+2,S_size(2)+2,3);
    temp_target(2:T_size(1)+1, 2:T_size(2)+1,:)=target;
    temp_mask(2:M_size(1)+1, 2:M_size(2)+1,:)=mask;
    temp_source(2:S_size(1)+1, 2:S_size(2)+1,:)=source;
    
    source=temp_source;
    mask=temp_mask;
    target=temp_target;
    Maskedge=1;
end
[mask_index_i,mask_index_j]=find(mask(:,:,1)==1);
[not_mask_index_i,not_mask_index_j]=find(mask(:,:,1)==0);

T_size=size(target);
M_size=size(mask);
S_size=size(source);
A_size=T_size(1)*T_size(2);

mask_ind=sub2ind([T_size(1),T_size(2)],mask_index_i,mask_index_j);
Msize=size(mask_ind);
M=Msize(1);

not_mask_ind=sub2ind(T_size,not_mask_index_i,not_mask_index_j);
NMsize=size(not_mask_ind);
NM=NMsize(1);

up_neighb_ind=sub2ind(T_size,mask_index_i-1,mask_index_j);
down_neighb_ind=sub2ind(T_size,mask_index_i+1,mask_index_j);
left_neighb_ind=sub2ind(T_size,mask_index_i,mask_index_j-1);
right_neighb_ind=sub2ind(T_size,mask_index_i,mask_index_j+1);

% x         y           v
% notmask, notmask -> 1
% mask, mask -> 4
% mask, neighbours(mask) -> -1
% O.W       ->          0

X=cat(1,not_mask_ind,mask_ind, mask_ind,mask_ind,mask_ind,mask_ind);% last four mask index for neghbours
Y=cat(1,not_mask_ind,mask_ind, up_neighb_ind,down_neighb_ind,left_neighb_ind,right_neighb_ind);
Values=cat(1,ones(NM,1), 4*ones(M,1), -1*ones(M,1), -1*ones(M,1), -1*ones(M,1), -1*ones(M,1));
A=sparse(X,Y,Values,A_size,A_size);

%b for 3 channels;

weight=1;

target1=target(:,:,1);
b1=target1(:);
target1=target1(:);
length=size(target1);
source1=source(:,:,1);
source1=source1(:);
b1(mask_ind)=4*source1(mask_ind)-source1(up_neighb_ind)-source1(down_neighb_ind)-source1(left_neighb_ind)-source1(right_neighb_ind);
b1_mixed=zeros(length(1),1);
b1_mixed(mask_ind)=4*target1(mask_ind)-target1(up_neighb_ind)-target1(down_neighb_ind)-target1(left_neighb_ind)-target1(right_neighb_ind);


%Comp=abs(cat(2,b1,b1_mixed)) ;
%[~,argmax]=max( Comp,[],2);
%b1(argmax(mask_ind)==2)=b1_mixed(argmax(mask_ind)==2);



b1(mask_ind)=(weight*b1(mask_ind)+(1-weight)*b1_mixed(mask_ind));

target2=target(:,:,2);
b2=target2(:);
target2=target2(:);
source2=source(:,:,2);
source2=source2(:);
b2(mask_ind)=4*source2(mask_ind)-source2(up_neighb_ind)-source2(down_neighb_ind)-source2(left_neighb_ind)-source2(right_neighb_ind);
b2_mixed=zeros(length(1),1);
b2_mixed(mask_ind)=4*target2(mask_ind)-target2(up_neighb_ind)-target2(down_neighb_ind)-target2(left_neighb_ind)-target2(right_neighb_ind);

%Comp=abs(cat(2,b2,b2_mixed)) ;
%[~,argmax]=max( Comp ,[],2);
%b2(argmax(mask_ind)==2)=b2_mixed(argmax(mask_ind)==2);
b2(mask_ind)=(weight*b2(mask_ind)+(1-weight)*b2_mixed(mask_ind));



target3=target(:,:,3);
b3=target3(:);
target3=target3(:);
source3=source(:,:,3);
source3=source3(:);
b3(mask_ind)=4*source3(mask_ind)-source3(up_neighb_ind)-source3(down_neighb_ind)-source3(left_neighb_ind)-source3(right_neighb_ind);
b3_mixed=zeros(length(1),1);
b3_mixed(mask_ind)=4*target3(mask_ind)-target3(up_neighb_ind)-target3(down_neighb_ind)-target3(left_neighb_ind)-target3(right_neighb_ind);

%Comp=abs(cat(2,b3,b3_mixed)) ;
%[~,argmax]=max( Comp ,[],2);
%b3(argmax(mask_ind)==2)=b3_mixed(argmax(mask_ind)==2);
b3(mask_ind)=(weight*b3(mask_ind)+(1-weight)*b3_mixed(mask_ind));





% Solving the system
x1=A\b1;
x2=A\b2;
x3=A\b3;

out1=reshape(x1,T_size(1),T_size(2));
out2=reshape(x2,T_size(1),T_size(2));
out3=reshape(x3,T_size(1),T_size(2));
out_target=zeros(T_size(1),T_size(2),3);
if(Maskedge==0)
    out_target(:,:,1)=out1;
    out_target(:,:,2)=out2;
    out_target(:,:,3)=out3;
else
    out_target=zeros(T_size(1)-2,T_size(2)-2,3);
    out_target(:,:,1)=out1(2:T_size(1)-1,2:T_size(2)-1,:);
    out_target(:,:,2)=out2(2:T_size(1)-1,2:T_size(2)-1,:);
    out_target(:,:,3)=out3(2:T_size(1)-1,2:T_size(2)-1,:);
end
output=out_target;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As explained on the web page, we solve for output by setting up a large
% system of equations, in matrix form, which specifies the desired value or
% gradient or Laplacian (e.g.
% http://en.wikipedia.org/wiki/Discrete_Laplace_operator)

% The comments here will walk you through a conceptually simple way to set
% up the image blending, although it is not necessarily the most efficient
% formulation. 

% We will set up a system of equations A * x = b, where A has as many rows
% and columns as there are pixels in our images. Thus, a 300x200 image will
% lead to A being 60000 x 60000. 'x' is our output image (a single color
% channel of it) stretched out as a vector. 'b' contains two types of known 
% values:
%  (1) For rows of A which correspond to pixels that are not under the
%      mask, b will simply contain the already known value from 'target' 
%      and the row of A will be a row of an identity matrix. Basically, 
%      this is our system of equations saying "do nothing for the pixels we 
%      already know".
%  (2) For rows of A which correspond to pixels under the mask, we will
%      specify that the gradient (actually the discrete Laplacian) in the
%      output should equal the gradient in 'source', according to the final
%      equation in the webpage:
%         4*x(i,j) - x(i-1, j) - x(i+1, j) - x(i, j-1) - x(i, j+1) = 
%         4*s(i,j) - s(i-1, j) - s(i+1, j) - s(i, j-1) - s(i, j+1)
%      The right hand side are measurements from the source image. The left
%      hand side relates different (mostly) unknown pixels in the output
%      image. At a high level, for these rows in our system of equations we
%      are saying "For this pixel, I don't know its value, but I know that
%      its value relative to its neighbors should be the same as it was in
%      the source image".

% commands you may find useful: 
%   speye - With the simplest formulation, most rows of 'A' will be the
%      same as an identity matrix. So one strategy is to start with a
%      sparse identity matrix from speye and then add the necessary
%      values. This will be somewhat slow.
%   sparse - if you want your code to run quickly, compute the values and
%      indices for the non-zero entries in A and then construct 'A' with a
%      single call to 'sparse'.
%      Matlab documentation on what's going on under the hood with a sparse
%      matrix: www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf  %% S(i(k),j(k)) = v(k)


%   reshape - convert x back to an image with a single call.
%   sub2ind and ind2sub - how to find correspondence between rows of A and
%      pixels in the image. It's faster if you simply do the conversion
%      yourself, though.
%   see also find, sort, diff, cat, and spy


