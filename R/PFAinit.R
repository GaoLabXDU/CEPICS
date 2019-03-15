#' Initialize a MATLAB server
#'
#' Initialize the MATLAB server on the local computer, add the path as MATLAB workspace, and create PFA source file for MATLAB on the path.
#'
#'You must have installed MATLAB on your computer.
#'
#' @param path A string, which means the work dictionary of MATLAB server. Default is the current R working directory getting from getwd().
#' @return If open MATLAB server successfully, return the access of MATLAB, the variable will be used as the parameter of runPFA.When you have done all of your work, you must close the access by \code{close(matlab)}.
#' @examples
#' \donttest{
#' matlab=PFAinit(path=getwd())
#' #now, we can run PFA on MATLAB server
#'
#' #TODO
#' #for example:
#' res=runPFA(data=list(COAD_Methy, COAD_miRNA, COAD_mRNA),maxk=10,matlab=matlab)
#'
#' #after working on MATLAB, we should close the MATLAB server.
#' close(matlab)
#' }
#'
#' @seealso
#' \code{\link{runPFA}}
#'
#' @export
PFAinit <- function(path=getwd()){
  options(warn =-1)
  MainPFA="
function Y=MainPFA(data)
  N=length(data);
  sample_num = size(data{1},2);


  d=[];
  eig_v={};
  u={};
  x={};
  for i=1:N
  [ui,eig_vi]=Algorithm1(data{i},sample_num);

  for di=1:sample_num
  if sum(eig_vi(1:di))/sum(eig_vi)>0.8
  break;
  end
  end
  d(i)=di;

  [ui,eig_vi]=Algorithm1(data{i},di);
  u{i}=ui;
  eig_v{i}=eig_vi;
  xi = ui'*data{i}; %% the local sample-spectrum
    x{i}=xi;
end

% capture the global sample-spectrum according to Algorithm_4
d_num = min(d);

iter_num =1000;
lam_1 = 1;

Y = Algorithm4(x, sample_num, iter_num,  lam_1, d_num,N);%% Y is the global sample-spectrum

Y=Y';
  end
"
A1="
function [disc_set,disc_value]=Algorithm1(Train_SET,Eigen_NUM)
 %% capture local sample-spectrum based on SVD

[NN,Train_NUM]=size(Train_SET);

if NN<=Train_NUM



R=Train_SET*Train_SET'/(Train_NUM-1);

[V,S]=Find_K_Max_Eigen1(R,Eigen_NUM);
disc_value=S;
disc_set=V;

else % small samples???svd


R=(Train_SET)'*Train_SET/(Train_NUM-1);

[V,S]=Find_K_Max_Eigen(R,Eigen_NUM);
disc_value=S;
disc_set=zeros(NN,Eigen_NUM);

Train_SET=Train_SET/sqrt(Train_NUM-1);
for k=1:Eigen_NUM
disc_set(:,k)=(1/sqrt(disc_value(k)))*Train_SET*V(:,k);
end

end

function [Eigen_Vector,Eigen_Value]=Find_K_Max_Eigen(Matrix,Eigen_NUM)

[NN,NN]=size(Matrix);
[V,S]=eig(Matrix); %Note this is equivalent to; [V,S]=eig(St,SL); also equivalent to [V,S]=eig(Sn,St); %

S=diag(S);
[S,index]=sort(S);

Eigen_Vector=zeros(NN,Eigen_NUM);
Eigen_Value=zeros(1,Eigen_NUM);

p=NN;
for t=1:Eigen_NUM
Eigen_Vector(:,t)=V(:,index(p));
Eigen_Value(t)=S(p);
p=p-1;
end

function [Eigen_Vector,Eigen_Value]=Find_K_Max_Eigen1(Matrix,Eigen_NUM)

[NN,NN]=size(Matrix);
[V,S]=eig(Matrix); %Note this is equivalent to; [V,S]=eig(St,SL); also equivalent to [V,S]=eig(Sn,St); %

S=diag(S);
[S,index]=sort(S);

Eigen_Vector=zeros(NN,NN);
Eigen_Value=zeros(1,NN);

p=NN;
for t=1:NN
Eigen_Vector(:,t)=V(:,index(p));
Eigen_Value(t)=S(p);
p=p-1;
end
"
A2="

function w_1= Algorithm2( sample_num, w, E, lam_1,N)
%% fix Y, solve W based Algorithm 2
for i=1:N
    E_or(i,:)=E{i};
end

E_or=reshape(E_or',N*sample_num  ,1);

E = (E_or)/sum(E_or);

[E_add,index]=sort(E);

p=inf;

for i=N*sample_num:-1:1
    o=(2*lam_1+sum(E_add(1:i,1)))/i - E_add(i);

    if o>=0
        p=i;
        break
    end
end

o=(2*lam_1+sum(E_add(1:p,1)))/p;
w(1:p,1)=(o-E_add(1:p,1))/(2*lam_1);
w((p+1):(N*sample_num))=0;

w_1 = -100*zeros(size(w,1),1);

w_1(index)=w;

end

"
A4="
 function Y = Algorithm4(x, sample_num, iter_num, lam_1, d_num,N)
    % x is the local sample-spectrum
% Y is the global sample-spectrum
wpre=ones(N*sample_num,1)./(N*sample_num);

Y_final = [];
w_final = [];

final_err= inf;
w={};
s={};
L={};
for iter=1:iter_num
M=0;
for i=1:N
a=(i-1)*sample_num+1;
b=i*sample_num;
w{i} = diag(sqrt(wpre(a:b,1)));
s{i} = sum(diag(((eye(sample_num)-(x{i})\\(x{i})))*((eye(sample_num)-(x{i})\\(x{i})))'));

            M=M+(1/s{i})*(w{i}*(eye(sample_num)-(x{i}*w{i})\\(x{i}*w{i})))*(w{i}*(eye(sample_num)-(x{i}*w{i})\\(x{i}*w{i})))';
end
[Y,~]=FindKMinEigen(M,d_num+1);
Y=Y(:,2:(d_num+1))';

        for i=1:N
            L{i} = (Y*w{i})/(x{i}*w{i});
        end

        err = sum( diag(Y*M*Y') )
if err-final_err<0
final_err = err;
Y_final = Y;
w_final = wpre;
else
break
end

E = computeerr( wpre, sample_num, Y,L, x, s ,N);

w{1}= Algorithm2( sample_num, wpre, E, lam_1,N);
wpre=w{1};

end

Y=Y_final;
wpre=w_final;
end
"
cp="
function E = compute_err( w, sample_num, Y,L, x, s ,N)

for j=1:N
E{j}=[];
end
for i=1:sample_num
for j=1:N

E{j} =[E{j}; (1/s{j})*w(sample_num*(j-1)+i,1)*(Y(:,i)-L{j}*x{j}(:,i))'*(Y(:,i)-L{j}*x{j}(:,i))];
end
end

end


"
fd="
function [Eigen_Vector,Eigen_Value]=Find_K_Min_Eigen(Matrix, Eigen_NUM)

[NN,NN]=size(Matrix);
[V,S]=eig(Matrix); %Note this is equivalent to; [V,S]=eig(St,SL); also equivalent to [V,S]=eig(Sn,St); %

S=diag(S);
[S,index]=sort(S);

Eigen_Vector=zeros(NN,Eigen_NUM);
Eigen_Value=zeros(1,Eigen_NUM);

p=1;
for t=1:Eigen_NUM
    Eigen_Vector(:,t)=V(:,index(p));
    Eigen_Value(t)=S(p);
    p=p+1;
end
"

cat(A1,file = "Algorithm1.m",sep='\n',append = FALSE)
cat(A2,file = "Algorithm2.m",sep='\n',append = FALSE)
cat(A4,file = "Algorithm4.m",sep='\n',append = FALSE)
cat(cp,file = "computeerr.m",sep='\n',append = FALSE)
cat(fd,file = "FindKMinEigen.m",sep='\n',append = FALSE)
cat(MainPFA,file = "MainPFA.m",sep='\n',append = FALSE)

  Matlab$startServer()
  matlab1 <- Matlab(port=9999)
  isopen<-open(matlab1)
  if(isopen==FALSE)
    return(NULL)
  evaluate(matlab1, paste("userpath('",path,"')",sep=""))
  return(matlab1)
}

