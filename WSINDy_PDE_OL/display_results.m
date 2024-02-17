figind=0;

%% plot basis fcn

if toggle_plot_basis_fcn
%     figind = figind +1;
%     figure(figind); 
    [x1,t1]=meshgrid(xs_obs{1}(1:2*ms(1,1)+1),xs_obs{end}(1:2*ms(1,2)+1));
    mesh(x1,t1,Cfs_t(1,:)'*Cfs_x(3,:))
    xlabel('x'); ylabel('t')
    drawnow
end

%% plot coefficients

if and(toggle_plot_weights,mod(tOL,toggle_plot_weights)==0)
%     figind = figind+1;
%     figure(figind);
    if isequal(class(axi),'double')
        true_weights = repmat(axi(axi~=0),1,tOL);
        learned_weights = cell2mat(cellfun(@(x) x(axi~=0),Ws,'uni',0));
        legs = tags_pde_G(and(axi==0,W~=0),:);
        learned_weights_FP = W(and(axi==0,W~=0),:);
        legs_FP = categorical(tags_pde_G(and(axi~=0,W~=0)));
    elseif isequal(class(axi),'cell')
        true_weights = cell2mat(cellfun(@(x) x(x~=0), axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));
        learned_weights = cell2mat(cellfun(@(x,y) x(y~=0), Ws, axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));
        legs = tags_pde_G(axi{1}~=0);
        learned_weights_FP = W(and(axi{(Kmem+1)/2+tOL-1}==0,W~=0),:);
        legs_FP = categorical(tags_pde_G(and(axi{(Kmem+1)/2+tOL-1}==0,W~=0)));
    end
    numterms = size(true_weights,1);
    numplots = numterms + 1;%~isempty(legs_FP);
    for j=1:numterms
        subplot(numplots,1,j)
        plot(1:tOL,true_weights(j,:), 'r.-', 1:tOL,learned_weights(j,:), 'b.--')
        xlim([1 min(2*tOL,length(xs_obs{end}))])
        legend({[legs{j}, ' true'],[legs{j}, ' learned']},'location','eastoutside')
    end
    subplot(numplots,1,numplots)
    bar(legs_FP,learned_weights_FP)
end

%% plot data

for jj=1:min(toggle_plot_sol,n)
%     figind = figind +1;
%     figure(figind); 
    colormap(turbo(50))
    foo=[repmat({':'},1,dim-1),{Kmem}];
    U_obs_gap=cellfun(@(x)x(foo{1:dim-1},1+tOL:Kmem+tOL),U_obs,'uni',0);
    if dim==3
        surf(xs_OL{1},xs_OL{2},U_obs_gap{jj}(foo{1:dim})', 'EdgeColor','none')
        view([0 90])       
        zlabel('$u$','interpreter','latex','fontsize',14)
        xlim([xs_OL{1}(1) xs_OL{1}(end)])
        ylim([xs_OL{2}(1) xs_OL{2}(end)])
        colorbar
    elseif and(dim==2,length(xs_obs{1})>1)
        plot(xs_OL{1},U_obs_gap{jj}(foo{1:dim})')
        xlim([xs_OL{1}(1) xs_OL{1}(end)])
    else
        if n==2
            plot(U_obs_gap{1},U_obs_gap{2},'.')
            xlim([min(U_obs{1}) max(U_obs{1})])
            ylim([min(U_obs{2}) max(U_obs{2})])
        elseif n==3
            plot3(U_obs_gap{1},U_obs_gap{2},U_obs_gap{3},'.')
            xlim([min(U_obs{1}) max(U_obs{1})])
            ylim([min(U_obs{2}) max(U_obs{2})])
            zlim([min(U_obs{3}) max(U_obs{3})])
        end
    end
    set(gca, 'TickLabelInterpreter','latex','fontsize',14)
    xlabel('$x$','interpreter','latex','fontsize',14)
    ylabel('$t$','interpreter','latex','fontsize',14)
end


%% plot data fft

if toggle_plot_fft>0
%     figind = figind +1;
%     figure(figind); 
    coords=1:dim;
    for j=1:dim-1
        coords=[coords;circshift(coords(end,:),-1)];
    end
    for j=1:dim
        coordsj = coords(j,:);
        dd=coordsj(1);
        U_obs_gap=cellfun(@(x)x(space_inds{:},tOL:tOL+Kmem-1),U_obs,'uni',0);
        Ufft = abs(fft(permute(U_obs_gap{1},coordsj)));
        Ufft = reshape(Ufft,size(Ufft,1),[]);
        Ufft = mean(Ufft(floor(end/2):end,:),2);
        L = length(Ufft)-1;
        ks = -L:L;
        Ufft = [Ufft; flipud(Ufft(1:end-1))]/max(Ufft);
        subplot(dim,1,j)
            semilogy(ks,Ufft)
            hold on
            if dd<dim
                Cfs_ffts_vec = fft([zeros(1,length(xs_OL{dd})-2*ms(1,1)-1) Cfs_x(1,:)]);
            else
                Cfs_ffts_vec = fft([zeros(1,length(xs_OL{dd})-2*ms(1,2)-1) Cfs_t(1,:)]);
            end
            Cfs_ffts_vec=abs(Cfs_ffts_vec(floor(end/2):end));
            Cfs_ffts_vec=[Cfs_ffts_vec fliplr(Cfs_ffts_vec(1:end-1))];
            Cfs_ffts_vec=Cfs_ffts_vec/max(Cfs_ffts_vec);
            semilogy(ks,Cfs_ffts_vec)
            hold off     
            ylim([min(Ufft)*0.1 max(Ufft)])
            legend({'$\mathcal{F}(U_d)$','$\mathcal{F}(\phi_d)$'},'interpreter','latex','fontsize',14)
            title(['coord',num2str(dd)])
    end
    xlabel('k')
end

drawnow
