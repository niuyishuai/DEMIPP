function setupfig(mode)
    % 设置坐标轴
    axes1 = gca; %axes('Parent',h);
    box(axes1,'on');
    set(axes1,'looseInset',[0 0 0 0]);
    
    % 设置其余坐标区属性
    set(axes1,'FontSize',16,'LineWidth',1.5);
    
    % labels
    xlabel('$v_1$','Interpreter','latex','FontWeight','bold','FontSize',23);
    ylabel('$v_2$','Interpreter','latex','FontWeight','bold','FontSize',23);
    
    if nargin<1
        mode=1;
    end
    switch mode
        case 1
            grid on;
            h=gcf;
            % 设置 figure属性
            set(h,'InvertHardcopy','off','PaperUnits','points',...
                'Color',[1 1 1],...
                'Renderer','painters',...
                'position',[100 300 800 530]);            
        case 2
            box off;
            zlabel('$\Pi$','Interpreter','latex','FontWeight','bold','FontSize',23);
            h=gcf;
            % 设置 figure属性
            set(h,'InvertHardcopy','off','PaperUnits','points',...
                'Color',[1 1 1],...
                'Renderer','painters',...
                'position',[1000 300 800 530]);
            view(3);
    end
end
