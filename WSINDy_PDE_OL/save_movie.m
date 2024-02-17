name = sprintf(['~/Desktop/wave2D_coeff_',date]);                          % Name for the video file
vidObj = VideoWriter(name);                                      % Creating video file
vidObj.FrameRate = 30;  open(vidObj);                            % Opening video file

fig202 = figure(202);                           		 % Initial plots
set(gcf,'position',[0 0 1000 800])

writeVideo(vidObj, getframe(fig202));
tic

legs = {'$\partial_{yy}u$','$\partial_{xx}u$','$u^3$'};
for tOL = 1:length(Ws)

    true_weights = cell2mat(cellfun(@(x) x(x~=0), axi((Kmem+1)/2:end),'uni',0));
    if isequal(class(axi),'double')
        true_weights = repmat(axi(axi~=0),1,tOL);
        learned_weights = cell2mat(cellfun(@(x) x(axi~=0),Ws,'uni',0));
        learned_weights_FP = W(and(axi==0,W~=0),:);
    elseif isequal(class(axi),'cell')
        learned_weights = cell2mat(cellfun(@(x,y) x(y~=0), Ws(1:tOL), axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));
        learned_weights_FP = cell2mat(cellfun(@(x,y) length(find(x(y==0))), Ws(1:tOL), axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));
%         learned_weights_FP = cell2mat(cellfun(@(x,y) tpscore(x,y), Ws(1:tOL), axi((Kmem+1)/2:(Kmem+1)/2+tOL-1),'uni',0));
    end
    numterms = size(true_weights,1);
    numplots = numterms+1;%~isempty(legs_FP);

    j=1;
    ax1=subplot(numplots,2,j*2);
    plot(1:tOL,true_weights(j,1:tOL), 'r.-', 1:tOL,learned_weights(j,:), 'b.--','linewidth',2)
    xlim([1 min(2*tOL,length(xs_obs{end}))])
    ylim([-0.5 4.5])
    legend({[legs{j}, ' true'],[legs{j}, ' learned']},'location','ne','fontsize',12,'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    p = get(gca,'Position');
    p([1 3]) = [0.4647 0.5];
    set(gca,'Position',p);


    j=2;
    ax2=subplot(numplots,2,j*2);
    plot(1:tOL,true_weights(j,1:tOL), 'r.-', 1:tOL,learned_weights(j,:), 'b.--','linewidth',2)
    xlim([1 min(2*tOL,length(xs_obs{end}))])
    ylim([-0.5 4.5])
    legend({[legs{j}, ' true'],[legs{j}, ' learned']},'location','ne','fontsize',12,'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    p = get(gca,'Position');
    p([1 3]) = [0.4647 0.5];
    set(gca,'Position',p);


    j=3;
    ax3=subplot(numplots,2,j*2);
    plot(1:tOL,true_weights(j,1:tOL), 'r.-', 1:tOL,learned_weights(j,:), 'b.--','linewidth',2)
    xlim([1 min(2*tOL,length(xs_obs{end}))])
    ylim([-10 0])
    legend({[legs{j}, ' true'],[legs{j}, ' learned']},'location','se','fontsize',12,'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    p = get(gca,'Position');
    p([1 3]) = [0.4647 0.5];
    set(gca,'Position',p);  

    ax4=subplot(numplots,2,numplots*2);
    plot(1:tOL,learned_weights_FP,'linewidth',2)
    xlim([1 min(2*tOL,length(xs_obs{end}))])
    ylim([0 40])
    xlabel('iter','Interpreter','latex')
    legend('false pos.','location','ne','fontsize',12,'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    p = get(gca,'Position');
    p([1 3]) = [0.4647 0.5];
    set(gca,'Position',p);

    subplot(numplots,2,1:2:numplots*2);
    imagesc(U_obs{1}(:,:,(Kmem+1)/2+tOL-1)')
%     contourf(U_obs{1}(:,:,(Kmem+1)/2+tOL-1)',50)
    axis equal
    xlim([0 length(xs_obs{1})])
    ylim([0 length(xs_obs{2})])
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
    colormap(turbo(50))
    colorbar('TickLabelInterpreter','latex','fontsize',13); caxis([-5 5])
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    p = get(gca,'Position'); p(1)=0;
    set(gca,"Position", p)
    title('$\partial_{tt}u = c(t)\Delta u - u^3$', 'interpreter','latex')
% drawnow
    writeVideo(vidObj, getframe(fig202));
end
toc
close(vidObj);