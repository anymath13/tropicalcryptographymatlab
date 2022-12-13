function [eigval,x,criticalcircuit]= policyIterationnonzeros(A)
%input:
%A : matriks ukuran nxnxd 
%d=dimensi
%tokenmaks=d-1
%eigval : eigen mode ukuran nx1
%x : Vektor eigen ukuran nx1


    e=-inf;
    n=size(A,1)
    d=size(A,3)
    tokenmaks=d-1
    x=zeros(n,1)
    eigval=zeros(n,1)
    if d==1 
        A1=A
        A=maxpluszeros(n,n)
        A(:,:,2)=A1
    end
    d=size(A,3)
    tokenmaks=d-1
    %cek matriks A mempunyai paling sedikit satu elemen tidak nol pada tiap baris
    Aom=maxpluszeros(n,n)
    for i=1:d
        Aom=add(Aom,A(:,:,i))
    end
    
    
    %Jika A0 zeros
    if max(A(:,:,1))==-inf 
            tanda=2;
    else
        tanda=1;
   
    %cek tidak ada cycle di A0
    
    
    %Tahap 1
    % Memilih sebarang policy pi
    pi=maxpluszeros(n,n)
    %banyak token
    mu=zeros(n,1)
    %pemetaan pi
    temp_pi=zeros(n,1)
    
    for j=1:n
        cari=find(Aom(j,:)~=e)
        for i=tanda:d
            if A(j,cari(1),i)~=e 
                pi(j,cari(1))=A(j,cari(1),i)
                temp_pi(j)=cari(1)
                mu(j)=i-1
                break;
            end
        end
    end
 
 waktu=0;   
 waktu2=0;
 Iterasi=0;
while 1    
    Iterasi=Iterasi+1
    %Tahap 3
    sim_pi=pi
    cek=maxpluszeros(n,1)
    
    banyakSimpankritis=0
    lintasankritis=[]
    while 1
        %menentukan sebarang cycle
       
        psi=cicle(pi)
        
        psi=[psi,psi(1)]
        psi100=psi
        psi100(length(psi100))=[]
        cekpanjang=length(psi100)
        if banyakSimpankritis==0 
            lintasankritis=psi100
        else
            if size(lintasankritis,2)<cekpanjang 
                lintasankritis=[lintasankritis,maxpluszeros(size(lintasankritis,1),cekpanjang-size(lintasankritis,2))]
            elseif size(lintasankritis,2)>cekpanjang
                psi100=[psi100,maxpluszeros(1,-cekpanjang+size(lintasankritis,2))]
            end
            lintasankritis=[lintasankritis;psi100]
        end
        criticalcircuit=lintasankritis
        banyakSimpankritis=banyakSimpankritis+1
       
        %menentukan eigen mode bar
        eigvalbar=[0,0]
        %disp('eigvalbar');
        for i=1:length(psi)-1
            eigvalbar(1)=eigvalbar(1)+pi(psi(i+1),psi(i))
            eigvalbar(2)=eigvalbar(2)+mu(psi(i+1))
            %disp(eigvalbar);
        end
        eigvalbar=eigvalbar(1)/eigvalbar(2)
        
        x1=x
        %pilih sebarang transisi
        pilih=psi(1)
        %update eigen mode untuk transisi pilihan
        eigval(pilih)=eigvalbar
        x1(pilih)=x(pilih)
        pi(pilih,temp_pi(pilih))=-inf
        newdirected=pilih
        %disp('pi sebelum terhubung');
        %disp(pi);
        cek(pilih)=1
        %disp('new directed');
        %disp(newdirected);
        %mendefinisikan nilai2 eigen dan vektor eigen untuk transisi yang lain yang terhubung terhadap transisi yang dipilih
        while 1
            directed=newdirected
            newdirected=[]
            for i=1:length(directed)
                cari=find(pi(:,directed(i))~=-inf)
                
                for j=1:length(cari)
                    if cek(cari(j))==-inf 
                        eigval(cari(j))=eigvalbar;
                        x1(cari(j))=pi(cari(j),temp_pi(cari(j)))-mu(temp_pi(cari(j)))*eigvalbar+x1(temp_pi(cari(j)))
                        newdirected=[newdirected,cari(j)]
                        pi(cari(j),temp_pi(cari(j)))=-inf
                        cek(cari(j))=1
                        %disp('pi cari terhubung');
                        %disp(pi);
                    end
                end
            end
            if isempty(newdirected)
                break;
            end
        end
       
        
        x=x1
        if max(max(pi))==-inf 
            break;
        end
    end
    x
    eigval
%    
    J=[]
%    
    K=zeros(n,1);
    I=[];
    L=zeros(n,1);
    token_K=zeros(n,1);
    token_L=zeros(n,1);
    
    
    for j=1:n
        Ktemp=[]
        Ltemp=[]
        tokenKtemp=[]
        Pj=find(Aom(j,:)~=-inf)
        eigpilih=eigval(Pj)
        maks=max(eigpilih)
        if (maks-eigval(j))>0.00000001
            J=[J,j]
        end
        f=find((abs(eigpilih-maks))==0)
        Pj=Pj(f)
        %penentuan K dan tokennya
        
        makss=-inf
        hit=1;
        sst=0;
        st=0;
        for s=tanda:d
            for ss=1:length(Pj)
                if A(j,Pj(ss),s)~=-inf 
%                   
                    banding=A(j,Pj(ss),s)-(s-1)*eigval(Pj(ss))+x(Pj(ss))
                    if banding>makss 
%                           
                            makss=banding
                            sst=ss
                            st=s
                    end
                end
            end
        end
        K(j)=Pj(sst)
        token_K(j)=st-1
        L(j)=Pj(sst)
        token_L(j)=st-1
     
        
        %Penentuan I
        if (makss-x(j))>0.00000001
            I=[I,j];
        end
        
        
    
    end
    
    pi=sim_pi;
   
    if isempty(I) && isempty(J) 

        return;
        
    elseif ~isempty(J)
        for ss=1:length(J)
            pi(J(ss),temp_pi(J(ss)))=-inf;
            mu(J(ss))=token_K(J(ss));
            temp_pi(J(ss))=K(J(ss));
            pi(J(ss),K(J(ss)))=A(J(ss),K(J(ss)),token_K(J(ss))+1);
           
        end
    else
        for ss=1:length(I)
            pi(I(ss),temp_pi(I(ss)))=-inf;
            mu(I(ss))=token_L(I(ss));
            temp_pi(I(ss))=L(I(ss));
            pi(I(ss),L(I(ss)))=A(I(ss),L(I(ss)),token_L(I(ss))+1);
            
        end
        
    end
  
    
end


end