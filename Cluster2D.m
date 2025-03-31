
 
 function [MClusterm1,MCluster0,MClusterp1]= Cluster2D(MAT)
dime=length(MAT);
m1=-1;
z0=0;
p1=1;
M_m1=find(MAT'==m1);
M_0=find(MAT'==z0);
M_p1=find(MAT'==p1);
M_adj_m1=zeros(length(M_m1));
M_adj_0=zeros(length(M_0));
M_adj_p1=zeros(length(M_p1));
MClusterm1=zeros(1,dime*dime);
MCluster0=zeros(1,dime*dime);
MClusterp1=zeros(1,dime*dime);

for i=1:dime
    im=i-1;
    if(im==0)
        im=dime;
    end
    ip=i+1;
    if(ip==dime+1)
        ip=1;
    end
    
        for j=1:dime
            
            jm=j-1;
            if(jm==0)jm=dime;end
            jp=j+1;
            if(jp==dime+1)jp=1;end
            
            if(MAT(i,j)==-1)
                
                iim1=find(M_m1==dime*(i-1)+j);
               
                if(MAT(i,j)==MAT(im,j))    
                    jjm1=find(M_m1==dime*(im-1)+j);
                    M_adj_m1(iim1,jjm1)=1;
                end
                
                if(MAT(i,j)==MAT(ip,j))
                    jjm1=find(M_m1==dime*(ip-1)+j);
                    M_adj_m1(iim1,jjm1)=1;
                end
                
                if(MAT(i,j)==MAT(i,jm))
                    jjm1=find(M_m1==dime*(i-1)+jm);
                    M_adj_m1(iim1,jjm1)=1;
                end
                    
                if(MAT(i,j)==MAT(i,jp))
                    jjm1=find(M_m1==dime*(i-1)+jp);
                    M_adj_m1(iim1,jjm1)=1;
                end
            end
            
              if(MAT(i,j)==0)
                
                ii0=find(M_0==dime*(i-1)+j);
              
                if(MAT(i,j)==MAT(im,j))
                    jj0=find(M_0==dime*(im-1)+j);
                    M_adj_0(ii0,jj0)=1;end
                    
                if(MAT(i,j)==MAT(ip,j))
                    jj0=find(M_0==dime*(ip-1)+j);
                    M_adj_0(ii0,jj0)=1;end
                    
                if(MAT(i,j)==MAT(i,jm))
                     jj0=find(M_0==dime*(i-1)+jm);
                    M_adj_0(ii0,jj0)=1;end
                    
                if(MAT(i,j)==MAT(i,jp))
                    jj0=find(M_0==dime*(i-1)+jp);
                    M_adj_0(ii0,jj0)=1;end
              end
               
               if(MAT(i,j)==1)
                   
                iip1=find(M_p1==dime*(i-1)+j);
                
                if(MAT(i,j)==MAT(im,j))
                    jjp1=find(M_p1==dime*(im-1)+j);
                    M_adj_p1(iip1,jjp1)=1;
                end
                    
                if(MAT(i,j)==MAT(ip,j))
                    jjp1=find(M_p1==dime*(ip-1)+j);
                    M_adj_p1(iip1,jjp1)=1;
                end
                               
                if(MAT(i,j)==MAT(i,jm))
                    jjp1=find(M_p1==dime*(i-1)+jm);
                    M_adj_p1(iip1,jjp1)=1 ;
                end
                    
                    
                if(MAT(i,j)==MAT(i,jp))
                   jjp1=find(M_p1==dime*(i-1)+jp);
                    M_adj_p1(iip1,jjp1)=1 ;
                end
               end
    
    end
end

Gm1 = graph(M_adj_m1);
G0 = graph(M_adj_0);
Gp1 = graph(M_adj_p1);

[clu_bm1, clu_Sm1] = conncomp(Gm1);
[clu_b0, clu_S0]  = conncomp(G0);
[clu_bp1, clu_Sp1]  = conncomp(Gp1);

%Lm1=unique(clu_Sm1);
%L0=unique(clu_S0);
%Lp1=unique(clu_Sp1);



for i=1:length(clu_Sm1)
MClusterm1(clu_Sm1(i))=MClusterm1(clu_Sm1(i))+1;
end

for i=1:length(clu_S0)
MCluster0(clu_S0(i))=MCluster0(clu_S0(i))+1;
end

for i=1:length(clu_Sp1)
MClusterp1(clu_Sp1(i))=MClusterp1(clu_Sp1(i))+1;
end

end


