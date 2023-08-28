function [Df,radius] = compute_pFq_circle(f,c,d,nxy,h,sing_info)
% INPUT: f: function on which to apply the integral operator
% a,b : The last elements of a abd of b will provide the exponents in the integrand
% nxy : domain: number of nodes in left and right in the x and bottom and
% top in the y directions.
% h: interval in x and in y
% sing_info: position and order of the singularity
% OUTPUT:
% Df: Hypergeometric function
% radius : radius of the circle inside which we perform a Taylor approx in
% the central region
% New Parameters:
% -Radius of the central region will now depend on alpha and beta
% -Paths can be adapted -- each corner must be far from a singularity
num_pdng_nodes = 10;%Number of padding nodes
sing_pos = inf;
al = c(end)-1;
bet = d(end)-c(end)-1;

Xinteger = nxy(1)-num_pdng_nodes:nxy(2)+num_pdng_nodes;
Yinteger = nxy(1)-num_pdng_nodes:nxy(2)+num_pdng_nodes;

[Xinteger_r,Xinteger_i] = meshgrid(Xinteger,Yinteger(end:-1:1));

xr = h*Xinteger_r;
xi = h*Xinteger_i;
Zinteger = Xinteger_r+1i*Xinteger_i;
z = h*Zinteger;

Zinteger_mat = Zinteger(:);
Zmat = z(:);

fvals = f(z); %We apply the operator to this function
if isempty(sing_info) == 0
    if sing_info(2)<0
        sing_pos = sing_info(1);
        fvals(z==sing_pos)=0;
    end
end

za=0; %Base point

Ib_vals= (z-za).^al; Ib_vals(z==za)=0;

%Nodes inside domain
ind_z_list = find(real(Zinteger_mat)>=nxy(1) & real(Zinteger_mat)<=nxy(2) & ...
                  imag(Zinteger_mat)>=nxy(3) & imag(Zinteger_mat)<=nxy(4));

%Global Operator matrices
G1 = sparse(length(ind_z_list),length(Zmat)); %Standard Corrections
G2 = sparse(length(ind_z_list),length(Zmat)); %Standard Corrections
G = sparse(length(ind_z_list),length(Zmat)); %Standard Corrections
G_a = sparse(length(ind_z_list),length(Zmat)); %Singularity at za=0
G_b = sparse(length(ind_z_list),length(Zmat)); %Singularity at zb=z0
G_center = sparse(length(ind_z_list),length(Zmat)); %Correction at the center

%Initializing derivatives
Df = zeros(length(ind_z_list(:)),1);

%Problematic Region
num_h_prob = 11;% Will now depend on alpha and beta


%%Computes the endpoint correction stencils
ns = 20; %Number of nodes in the stencil
t=linspace(0,2*pi,ns+1)'; t(end)=[];

Ra = 4; %Radius of stencil around za
Rb = 4; %Radius of stencil around zb
Rc = 4; %Radius of stencil around corner
pa = 3; %Node at which the TR starts from za
pb = 3; %Node at which the TR starts from zb
pc = 3; %Node at which the TR starts from a corner


ind_a = find(ismember(Zinteger_mat,(round(Ra*cos(t))+1i*round(Ra*sin(t)))));
ind_b = find(ismember(Zinteger_mat,(round(Rb*cos(t))+1i*round(Rb*sin(t)))));
ind_c = find(ismember(Zinteger_mat,(round(Rc*cos(t))+1i*round(Rc*sin(t)))));

na = length(ind_a);
nb = length(ind_b);
nc = length(ind_c);

Zinteger_sa = Zinteger_mat(ind_a);
Zinteger_sb = Zinteger_mat(ind_b);
Zinteger_sc = Zinteger_mat(ind_c);

cw_lft = stencil_circle(0, 1, pc,Zinteger_sc); %Correction at Singularity-free ends
cw_rt  = stencil_circle(0,-1, pc,Zinteger_sc);
cw_top = stencil_circle(0,-1i,pc,Zinteger_sc);
cw_bot = stencil_circle(0, 1i,pc,Zinteger_sc);

cwa_lft = stencil_circle(al, 1, pa,Zinteger_sa); %Correction at origin
cwa_rt  = stencil_circle(al,-1, pa,Zinteger_sa);
cwa_top = stencil_circle(al,-1i,pa,Zinteger_sa);
cwa_bot = stencil_circle(al, 1i,pa,Zinteger_sa);

cwb_lft = stencil_circle(bet, 1, pb,Zinteger_sb); %Correction at evaluation point
cwb_rt  = stencil_circle(bet,-1, pb,Zinteger_sb);
cwb_top = stencil_circle(bet,-1i,pb,Zinteger_sb);
cwb_bot = stencil_circle(bet, 1i,pb,Zinteger_sb);


r = min(num_h_prob,sing_pos/h);%radius for central node algorithm
ind_center = find(ismember(Zinteger_mat,(round(r*cos(t))+1i*round(r*sin(t)))));
radius = r*h; %radius at which we transition methods


%Position of the origin in the matrix
[cntr_y,cntr_x]=find(Zinteger==0);

[row,col] = ind2sub(size(z),ind_z_list(:));

for ind_z=1:length(ind_z_list)

    Wtr = zeros(size(z)); %Weight matrix (Trapezoidal Rule)
    WaVec = zeros(size(Zmat));
    WbVec = zeros(size(Zmat));
    WcVec = zeros(size(Zmat));

    x0 = xr(row(ind_z),col(ind_z));%Evaluation point zb = x0 + i y0
    y0 = xi(row(ind_z),col(ind_z));
    zb=x0+1i*y0;

    if (abs(zb)<=radius)

        G_center(ind_z,ind_center)=central_weights(al,bet,zb,Zmat(ind_center));

    else

        ind_x0 = round(x0/h);ind_y0 = round(y0/h);%Eval point coordinates

        delta_Za = 10;%safe distance from Za
        delta_Zb = 10;%safe distance from Zb
        delta_Zc = 0;%safe distance from Zc (for cases s.a. f(z)=(z-1)^(-a))

        flagb=0; %if flagb == 1: will relocate the branch cut associated with Zb.
        flag_neg=0; %if flag_neg == 1: allow to integrate along neg real axis in the lower branch.
        if ind_x0==0 | ind_y0==0
            flagb = 0;
            path = [0 0; ind_x0 ind_y0];
        elseif ind_x0>=-delta_Za & abs(ind_y0)>delta_Zb
            flagb = 1;
            if ind_y0<0
                flag_neg=1;
            end
            path = [0 0; -delta_Za-delta_Zb 0;-delta_Za-delta_Zb ind_y0;ind_x0 ind_y0];
        elseif ind_x0>=0 & abs(ind_y0)<=delta_Zb
            flagb=1;
            path = [0 0; 0 sign(ind_y0)*(delta_Za+delta_Zb);ind_x0 sign(ind_y0)*(delta_Za+delta_Zb); ind_x0 ind_y0];
        elseif ind_x0<0 & ind_y0>delta_Zb %TEST POS
            flagb = 0;
            path = [0 0; ind_x0 0; ind_x0 ind_y0];
        elseif ind_x0<0 & ind_y0<-delta_Zb %TEST NEG
            flag_neg=1;
            flagb = 0;
            path = [0 0; ind_x0 0; ind_x0 ind_y0];
        elseif ind_x0<0 & abs(ind_y0)<=delta_Zb
            flagb = 0;
            path = [0 0; 0 sign(ind_y0)*(delta_Za+delta_Zb);ind_x0 sign(ind_y0)*(delta_Za+delta_Zb); ind_x0 ind_y0];
        end

        if flagb==1
            Ia_vals = (zb-z).^bet;
        else
            Ia_vals = (z-zb).^bet;
        end

        Ia_vals(z==zb)=0;
        Iab_vals = Ia_vals.*Ib_vals;

        for idx_pth = 2:size(path,1) %For each segment

            % Get TR weights
            stpx = sign(path(idx_pth,1)-path(idx_pth-1,1));dx = stpx+(stpx==0);
            stpy = sign(path(idx_pth,2)-path(idx_pth-1,2));dy = stpy+(stpy==0);
            path_yi=cntr_y-((path(idx_pth-1,2)):dy:(path(idx_pth,2)));
            path_xi=cntr_x+((path(idx_pth-1,1)):dx:(path(idx_pth,1)));

            dir = stpx+1i*stpy; dir = dir/abs(dir);
            h_dir = h*dir;

            corr_TR = ones(size(z(path_yi,path_xi)));

            if flag_neg==1 & abs(Zinteger(path_yi(1),path_xi(1)))<eps & abs(imag(Zinteger(path_yi(end),path_xi(end))))<eps
                corr_TR = corr_TR*exp(-1i*2*pi*al);
            end

            tr_path = h_dir*corr_TR;

            if idx_pth==2 & idx_pth == size(path,1)
                tr_path(1:pa)=0;
                tr_path(end+1-pb:end)=0;
            elseif idx_pth==2
                tr_path(1:pa)=0;
                tr_path(end+1-pc:end)=0;

            elseif idx_pth == size(path,1)
                tr_path(1:pc)=0;
                tr_path(end+1-pb:end)=0;
            else
                tr_path(1:pc)=0;
                tr_path(end+1-pc:end)=0;
            end

            Wtr(path_yi,path_xi)=Wtr(path_yi,path_xi)+tr_path;%TR weights

            % Get start correction weights
            if idx_pth==2 %start

                imag_zsa = h*imag(Zinteger_sa);

                sheet_corra = ones(na,1); %Correction across the branch cuts
                if y0<=2*Ra*h & y0>0 & x0>=0 & flagb==0
                    sheet_corra(imag_zsa>=y0)=exp(-2i*pi*bet);
                elseif y0>=(-2*Ra*h) & y0<=0 & x0>=0 & flagb==0
                    sheet_corra(imag_zsa<y0)=exp(2i*pi*bet);
                end

                if flag_neg==1
                    sheet_corra = sheet_corra*exp(-2i*pi*al);
                end


                WaVec(ind_a)=WaVec(ind_a)+(h_dir)^(al+1)*sheet_corra.*( ...
                    (stpy==0 & stpx>0)*cwa_lft+...%left start correction
                    (stpy==0 & stpx<0)*cwa_rt+...%lft end corr
                    (stpy>0 & stpx==0)*cwa_bot+...%bottom start correction
                    (stpy<0 & stpx==0)*cwa_top);
            else

                ind_cj = find(ismember(Zinteger_mat,Zinteger(path_yi(1),path_xi(1))+Zinteger_sc));
                Zc = Zinteger(path_yi(1),path_xi(1));

                sheet_corrc = ones(nc,1); %Correction across the branch cuts
                if imag(Zc)<=(2*Rc) & imag(Zc)>=0 & real(Zc)<0 & flag_neg==0
                    sheet_corrc(imag(Zinteger_mat(ind_cj))<0)=exp(2i*pi*al);
                elseif (imag(Zc)>=(-2*Rb) & imag(Zc)<0 & real(Zc)<0 & flag_neg==0) | (imag(Zc)==0 & real(Zc)<0 & flag_neg==1)
                    sheet_corrc(imag(Zinteger_mat(ind_cj))>=0)=exp(-2i*pi*al);
                end
                WcVec(ind_cj)=WcVec(ind_cj)+h_dir*sheet_corrc.*( ...
                    (stpy==0 & stpx>0)*cw_lft+...%left start correction
                    (stpy==0 & stpx<0)*cw_rt+...%right start correction
                    (stpy>0 & stpx==0)*cw_bot+...%bottom start correction
                    (stpy<0 & stpx==0)*cw_top);%top start correction

            end

            % Get end correction
            if idx_pth == size(path,1) %Last segment

                trnslted_Zinteger_sb = Zinteger(path_yi(end),path_xi(end))+Zinteger_sb;
                ind_b = find(ismember(Zinteger_mat,trnslted_Zinteger_sb));
                imag_zsb = imag(trnslted_Zinteger_sb);

                sheet_corrb = ones(nb,1); %Correction across the branch cuts
                if y0<=(Rb*h) & y0>=0 & x0<0
                    sheet_corrb(imag_zsb<0)=exp(2i*pi*al);
                elseif y0>=(-Rb*h) & y0<0 & x0<0
                    sheet_corrb(imag_zsb>=0)=exp(-2i*pi*al);
                    elseif y0>=(-2*Rb*h) & y0<=0 & x0>1 & length(c)==2
                    sheet_corrb(imag_zsb>0)=exp(2i*pi*c(1));
                    elseif y0<=(2*Rb*h) & y0>0 & x0>1 & length(c)==2
                    sheet_corrb(imag_zsb<=0)=exp(-2i*pi*c(1));
                end

                WbVec(ind_b)=WbVec(ind_b)+((-h_dir)^(bet+1))*sheet_corrb.*(...
                    (stpy==0 & stpx<0)*cwb_lft+...%right end correction
                    (stpy==0 & stpx>0)*cwb_rt+....%lft end corr top
                    (stpy<0 & stpx==0)*cwb_bot+....%top end correction
                    (stpy>0 & stpx==0)*cwb_top);%bottom end correction

            else

                %Find be an easier way since on a regular grid
                ind_cj = find(ismember(Zinteger_mat,Zinteger(path_yi(end),path_xi(end))+Zinteger_sc));
                Zc = Zinteger(path_yi(end),path_xi(end));

                sheet_corrc = ones(nc,1); %Correction across the branch cuts
                if imag(Zc)<=(2*Rc) & imag(Zc)>=0 & real(Zc)<0 & flag_neg==0
                    sheet_corrc(imag(Zinteger_mat(ind_cj))<0)=exp(2i*pi*al);
                elseif (imag(Zc)>=(-2*Rb) & imag(Zc)<0 & real(Zc)<0)  | (imag(Zc)==0 & real(Zc)<0 & flag_neg==1)
                    sheet_corrc(imag(Zinteger_mat(ind_cj))>=0)=exp(-2i*pi*al);
                end



                WcVec(ind_cj)=WcVec(ind_cj)+h_dir*sheet_corrc.*( ...
                    (stpy==0 & stpx<0)*cw_lft+... %left end correction
                    (stpy==0 & stpx>0)*cw_rt+... %right end correction
                    (stpy<0 & stpx==0)*cw_bot+...%bottom end correction
                    (stpy>0 & stpx==0)*cw_top); %top end correction
            end


        end



        % Transform the matrices W with Wb into the j^th row of global
        % matrix G
        [~,Gpos1,Gval1] = find(((WcVec).*Iab_vals(:)).');
        G1(ind_z_list(ind_z),Gpos1)=Gval1;
        [~,Gpos2,Gval2] = find(((Wtr(:)).*Iab_vals(:)).');
        G2(ind_z_list(ind_z),Gpos2)=Gval2;

        [~,Gpos,Gval] = find(((WcVec+Wtr(:)).*Iab_vals(:)).');
        G(ind_z_list(ind_z),Gpos)=Gval;
        [~,Gpos_b,Gval_b] = find((WbVec.*Ib_vals(:)).');
        G_b(ind_z_list(ind_z),Gpos_b)=Gval_b;
        [~,Gpos_a,Gval_a] = find((WaVec.*Ia_vals(:)).');
        G_a(ind_z_list(ind_z),Gpos_a)=Gval_a;

    end

    if angle(x0+1i*y0)>0 & flagb==0
        Cb(ind_z) = exp(1i*pi*bet);
        Ca(ind_z) = exp(1i*pi*bet);
    elseif  flagb==0
        Cb(ind_z) = exp(-1i*pi*bet);
        Ca(ind_z) = exp(-1i*pi*bet);
    end

    if flagb==1;
        if ind_y0<0 & ind_y0>=-10
            Cb(ind_z) = exp(1i*pi*bet);
        elseif ind_y0<0 & ind_y0<=-10
            Cb(ind_z) = exp(-1i*pi*bet);
        elseif ind_y0>0 & ind_y0<=10
            Cb(ind_z) = exp(-sign(y0)*1i*pi*bet);
        elseif ind_y0>0 & ind_y0>10
            Cb(ind_z) = exp(-sign(y0)*1i*pi*bet);
        end
        Ca(ind_z) = 1;
    end
end
Dmat = G_center+diag(Cb)*(G_b(ind_z_list,:))+diag(Ca)*G2(ind_z_list,:)-diag(Ca)*G1(ind_z_list,:)-diag(Ca)*G_a(ind_z_list,:);

% figure(10)
% [S1,S2]=size(Dmat);
% [Di,Dj,Dv] = find(Dmat);
% scatter(Dj,S1-Di,20,log10(abs(Dv)))

Int= Dmat*fvals(:);
Df = (gamfun(d(end))/(gamfun(c(end))*gamfun(d(end)-c(end))))*(Int).*(z(ind_z_list(:)).^(1-d(end)));
