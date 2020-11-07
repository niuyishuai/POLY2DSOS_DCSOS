function setupfig(pos)
    if nargin<1
        pos=1;
    end
    sz=get(0,'ScreenSize');
    H=sz(4); % screen height
    W=sz(3); % screen width
    C=[W/2,H/2]; % screen center
    figsize=[800,500]; % figure size
    marginsize=abs((C-figsize)/2); % margin size
    % [1600 1000 800 530]
    switch pos
        case 1 % left up
            position=[marginsize(1),C(2)-50,figsize];
        case 2 % right up
            position=[C(1)+marginsize(1),C(2)-2*marginsize(2),figsize];
        case 3 % left down
            position=[marginsize(1),2*marginsize(2),figsize];
        case 4 % right down
            position=[C(1)+marginsize(1),2*marginsize(2),figsize];
    end
    % set axis
    axes1 = gca;
    box(axes1,'on');
    set(axes1,'looseInset',[0 0 0 0]);
    
    % set axis fonts and lines
    set(axes1,'FontSize',16,'LineWidth',1.5);
    
    grid off;
    h=gcf;
    % set fig pos
    set(h,'InvertHardcopy','off','PaperUnits','points',...
        'Color',[1 1 1],...
        'Renderer','painters',...
        'position',position);
    
end
