
function plot_network(net_vect,sp_names,sp_cols,intDir,intThick,ttl)
    ns = length(sp_names);
    mdl = reshape(net_vect,[ns ns])';

    G = digraph(mdl);
    h = plot(G);
    
    load('colormaps.mat')
    colormap(myredblue);
    h.EdgeCData = intDir;
    h.LineWidth = intThick.*3;
    h.NodeColor = sp_cols;


    h.NodeLabel = sp_names;
    h.MarkerSize = 20;
    h.ArrowSize = 15;
    h.ArrowPosition = 0.9;
    h.NodeFontSize = 12;

    title(ttl)
end