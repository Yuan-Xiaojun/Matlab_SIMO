function  A = Constell_Modulate(B, qam)

[M,N] = size(B);
A = zeros(M/sqrt(qam), N);



%%
     if qam == 4
        for  n=1 : N
            for  m=1 : M/sqrt(qam)
                pack=[B((2*m-1),n) B(2*m,n)];
                if  pack==[1 1] 
                    A(m,n)=1+1i;
               elseif  pack==[0 1] 
                    A(m,n)=-1+1i;
               elseif  pack==[0 0] 
                    A(m,n)=-1-1i;       
                else
                    A(m,n)=1-1i;
                end
            end
        end
        A = A./sqrt(2);
     elseif qam == 16
             for  n=1 : N
                for  m=1 : M/sqrt(qam)
                    pack=[B((4*m-3),n) B((4*m-2),n) B((4*m-1),n) B(4*m,n)];
                    if  pack == [1 0 1 0] 
                        A(m,n) = 1+1i;
                   elseif  pack == [0 1 1 0] 
                        A(m,n) = -1+1i;
                   elseif  pack == [0 1 0 1] 
                        A(m,n) = -1-1i;
                    elseif  pack == [1 0 0 1] 
                        A(m,n) = 1-1i;
                   elseif  pack == [1 1 1 0] 
                        A(m,n) = 3+1i;
                   elseif  pack == [0 1 1 1] 
                        A(m,n) = -1+3i;
                   elseif  pack == [0 0 0 1] 
                        A(m,n) = -3-1i;
                   elseif  pack == [1 0 0 0] 
                        A(m,n) = 1-3i;
                   elseif  pack == [1 1 1 1] 
                        A(m,n) = 3+3i;
                   elseif  pack == [0 0 1 1] 
                        A(m,n)=-3+3i;
                   elseif  pack == [0 0 0 0] 
                        A(m,n) = -3-3i;
                   elseif  pack == [1 1 0 0] 
                        A(m,n) = 3-3i;
                   elseif  pack == [1 0 1 1] 
                        A(m,n) = 1+3i;
                   elseif  pack == [0 0 1 0] 
                        A(m,n) = -3+1i;
                   elseif  pack == [0 1 0 0] 
                        A(m,n) = -1-3i;
                     else   
                        A(m,n) = 3-1i;
                    end
                end  
             end 
             A = A./sqrt(10);  
      end
end