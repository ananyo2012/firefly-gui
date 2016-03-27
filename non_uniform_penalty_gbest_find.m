%d=0.5*[0.470 1.322 2.263 3.178 4.142 5.369 6.212 7.135 8.313 9.794 11.192 12.792 14.360 15.960];
%d=[0.224 .573 .973 1.375 1.788 2.161 2.607 3.041 3.632 4.125 4.661 5.458 6.203 7.011];
x=gbest_find;
%x=nds(8,:);
%x=diff(d);
FNBWD=8.35;
tol=(12.0/100)*FNBWD;
P1=10^6;
P2=10^6;
lb1=-90;ub1=90; 
 step=.1;ed=(ub1-lb1)/step;
 NULL(1)=30;NULL(2)=32.5;NULL(3)=35;
 NULL(1)=((NULL(1)-lb1)/step)+1;
 NULL(2)=((NULL(2)-lb1)/step)+1;
 NULL(3)=((NULL(3)-lb1)/step)+1;
 NULLD=-60;
  N=14; %Number of element
  d(1)=x(1);
  ii=1:N;s(ii)=1;%Uniform excitation current
  for ii=2:N
      d(ii)=d(ii-1)+x(ii);
  end
   jj=0;
       
    for theta=lb1:step:ub1
        %u1=u-eps;
        %ele=(cos(0.5*pi*u1))/((1-u1^2)^0.5);
        
          jj=jj+1; 
         for ii=1:N
              
             pencil(ii)=2*s(ii)*cos(2*pi*d(ii)*sin(theta*pi/180));
             
          end
        pencil_array(jj)=abs(sum(pencil));
       
    end
    npencil_array=pencil_array/max(pencil_array);
    
   
    db_npencil_array=20*log10(npencil_array);
    
      
   minVal  = sign(diff([inf npencil_array inf]));
 indMin  = find(diff(minVal+(minVal==0))==2);      % Indices FF nulls

 %--------------Find all FF peaks----------------------------------
 
 indPeaks = find(diff(minVal)<0);           % Indices FF peaks

 [peakLevel indP] = sort(npencil_array(indPeaks), 'descend');
 

  indMax     = indPeaks(indP(1));
 indNullL   = find (indMin<indMax,1, 'last');  % Index first null mainbeam lower

 indNullU   = find (indMin>indMax,1, 'first');  % Index first null mainbeam upper
 
 FNU=indMin(indNullU);
 FNL=indMin(indNullL);
%  FNU=943;
% FNL=859;
 FNBWU=(FNU-1)*step+lb1;
 FNBWL=lb1+(FNL-1)*step;
 FNBW=abs(FNBWU-FNBWL)
%  FNBW= 2*FNBWU;
 
 
 
 %-------------Find indices all SLL directions--------------------
 
 
 indSLL  = [1:FNL-1 FNU+1:ed+1];
 

 max_SLL=max(db_npencil_array(indSLL))%Maximim value SLL
% max_SLL  = 20*log10(peakLevel(2));    % Maximum value SLL  
%   if max_SLL<-25,err=0;else err=abs(max_SLL+25);end 
 if db_npencil_array(NULL(1))<NULLD, err1=0;else err1=abs(db_npencil_array(NULL(1))-NULLD);end
 if db_npencil_array(NULL(2))<NULLD, err2=0;else err2=abs(db_npencil_array(NULL(2))-NULLD);end
 if db_npencil_array(NULL(3))<NULLD, err3=0;else err3=abs(db_npencil_array(NULL(3))-NULLD);end
 err4=abs(FNBW-FNBWD);
% if FNBW>=FNBWD && FNBW<9.35, err4=0; else err4=abs(FNBW-9.35);end
  %r=max(0,db_npencil_array(NULL(1))-NULLD)+max(0,db_npencil_array(NULL(2))-NULLD)+max(0,db_npencil_array(NULL(3))-NULLD);
  NULL_depth_30=db_npencil_array(NULL(1))
  NULL_depth_32=db_npencil_array(NULL(2))
  NULL_depth_35=db_npencil_array(NULL(3))
  
  cost=max_SLL+P1*max(0,abs(FNBW-FNBWD)-tol);
  %cost=max_SLL+err1+err2+err3+err4;
%   if cost>0, 
%   cost=0; 
%   else
%       cost =cost;
%   end
distance_from_orgin=d
  figure(2);
  theta=lb1:step:ub1;
  plot(theta,db_npencil_array);
   