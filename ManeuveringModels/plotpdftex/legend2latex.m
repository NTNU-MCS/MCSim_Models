% LEGEND2LATEX(FIGURE_HANDLE) Converts the present legends in the figure to
% legends that can be postprocessed using LaPrint.
% 
% Note that all earlier legend2latex objects are removed.

function out = legend2latex(figure_handle)
% gets all legend entries and converts these to a cell
leghs = findobj(figure_handle,'Type','Axes','Tag','legend');
if isempty(leghs)
    error('No legend found! Are you sure you selected the figure (NOT the axes), i.e. gcf?')
elseif length(leghs) > 1
    warning('More than one legend found!')  
end

for j = 1:length(leghs)
    legh = leghs(j);
    
    % defaults
    object_opts = '''VerticalAlignment'', ''top'', ''HorizontalAlignment'', ''right''';
    legend_percentual_x = 0.98;
    legend_percentual_y = 0.98;
    
    % swithc off interpreter and get the legend variables
    set(legh, 'interpreter', 'none');
    legendstring = get(legh,'String');
    location = get(legh,'Location');
    
    % translate the location into a relative cordinate
    corner_dist = 0.05;
    switch lower(location)
        case 'north'
            legend_percentual_x = 0.5;
            legend_percentual_y = 1-corner_dist;
            object_opts = '''VerticalAlignment'', ''top'', ''HorizontalAlignment'', ''center''';
        case 'south'
            legend_percentual_x = 0.5;
            legend_percentual_y = corner_dist;
            object_opts = '''VerticalAlignment'', ''bottom'', ''HorizontalAlignment'', ''center''';
        case 'southeast'
            legend_percentual_x = 1-corner_dist;
            legend_percentual_y = corner_dist;
            object_opts = '''VerticalAlignment'', ''bottom'', ''HorizontalAlignment'', ''right''';
        case 'southwest'
            legend_percentual_x = corner_dist;
            legend_percentual_y = corner_dist;
            object_opts = '''VerticalAlignment'', ''bottom'', ''HorizontalAlignment'', ''left''';
        case 'northeast'
            legend_percentual_x = 1-corner_dist;
            legend_percentual_y = 1-corner_dist;
            object_opts = '''VerticalAlignment'', ''top'', ''HorizontalAlignment'', ''right''';
        case 'northwest'
            legend_percentual_x = corner_dist;
            legend_percentual_y = 1-corner_dist;
            object_opts = '''VerticalAlignment'', ''top'', ''HorizontalAlignment'', ''left''';
        case 'east'
            legend_percentual_x = 1-corner_dist;
            legend_percentual_y = 0.5;
            object_opts = '''VerticalAlignment'', ''middle'', ''HorizontalAlignment'', ''right''';
        case 'west'
            legend_percentual_x = corner_dist;
            legend_percentual_y = 0.5;
            object_opts = '''VerticalAlignment'', ''middle'', ''HorizontalAlignment'', ''left''';
        otherwise
            warning('legend value unknown')
    end
    
    % get the corresponding lines and colours
    ud = get(legh(1),'UserData');
    plots = get(ud.handles(1),'Parent');
    linecolors = zeros(length(ud.handles), 3);
    linewidthn = zeros(length(ud.handles), 1);
    linewidths = cell(1,length(ud.handles));
    linestyles = cell(1,length(ud.handles));
    markers = cell(1,length(ud.handles));
    markers_edge = cell(1,length(ud.handles));
    MarkerEdgeColor = cell(1,length(ud.handles));
    MarkerFaceColor = cell(1,length(ud.handles));
    noedge = zeros(length(ud.handles), 1);
    patch = zeros(length(ud.handles), 1);
    patchColor = zeros(length(ud.handles), 3);
    patchEdge = zeros(length(ud.handles), 3);
    is_stem = 0;
    for i = 1:length(ud.handles)
    if strcmp(get(ud.handles(i), 'Type'),'hggroup')
        try get(ud.handles(i), 'FaceColor')
            catch Merr
                Merr;
        end
        if ~exist('Merr','var')
        patch(i) = 1;
        patchColor(i,:) = get(ud.handles(i), 'FaceColor');
        if isnumeric(patchColor(i,:)) && length(patchColor(i,:)) == 3 && max(patchColor(i,:)) <= 1
            patchColor(i,:) = patchColor(i,:);
        else
            warning('FaceColor value not recognized')
        end
        patchEdge(i,:) = get(ud.handles(i), 'EdgeColor');
        if isnumeric(patchEdge(i,:)) && length(patchEdge(i,:)) == 3 && max(patchEdge(i,:)) <= 1
            patchEdge(i,:) = patchEdge(i,:);
        else
            warning('EdgeColor value not recognized')
        end
        else
            is_stem = 1;
        end
    end
    if strcmp(get(ud.handles(i), 'Type'),'line') || is_stem
        linecolors(i,:) = get(ud.handles(i), 'Color');
        if isnumeric(linecolors(i,:)) && length(linecolors(i,:)) == 3 && max(linecolors(i,:)) <= 1
            linecolors(i,:) = linecolors(i,:);
        else
            warning('LineColor value not recognized')
        end
        linewidthn(i,:) = get(ud.handles(i), 'LineWidth');
        if isnumeric(linewidthn(i,:)) && length(linewidthn(i,:)) == 1
            linewidths{i} = num2str(linewidthn(i,:));
        else
            warning('LineWidth value not recognized')
        end
        linestyles{i} = get(ud.handles(i), 'LineStyle');
        switch linestyles{i}
            case '-'
               linestyles{i} = ['\rule[0.5ex]{3em}{',linewidths{i},'pt}'];
            case '--'
                linestyles{i} = ['\rule[0.5ex]{0.625em}{',linewidths{i},'pt}\hspace{0.55em}\rule[0.5ex]{0.625em}{',linewidths{i},'pt}\hspace{0.55em}\rule[0.5ex]{0.625em}{',linewidths{i},'pt}'];
            case ':'
                linestyles{i} = ['\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}\hspace{0.2em}\rule[0.5ex]{0.2em}{',linewidths{i},'pt}'];
            case '-.'
                linestyles{i} = ['\rule[0.5ex]{0.625em}{',linewidths{i},'pt}\hspace{0.16em}\rule[0.5ex]{0.25em}{',linewidths{i},'pt}\hspace{0.16em}\rule[0.5ex]{0.625em}{',linewidths{i},'pt}\hspace{0.16em}\rule[0.5ex]{0.25em}{',linewidths{i},'pt}\hspace{0.16em}\rule[0.5ex]{0.625em}{',linewidths{i},'pt}'];
            case 'none'
                linestyles{i} = '';
            otherwise
                warning('Linestyle value not recognized');
                linestyles{i} = ['\rule[0.5ex]{3em}{',linewidths{i},'pt}'];
        end
        markers{i} = get(ud.handles(i), 'Marker');
        switch markers{i}
            case 'none'
                markers{i} = '';
            case '+'
                markers{i} = '$\mathbf{+}$';
            case 'o'
                markers{i} = '\put(-0.5,-3){\huge{$\bullet$}}';
                markers_edge{i} = '\put(-0.5,-3){\huge{$\circ$}}';
            case '*'
                markers{i} = '\put(0,-7.5){\huge{*}}';
            case '.'
                markers{i} = '\put(1,0){\small{$\bullet$}}';
            case 'x'
                markers{i} = '$\mathbf{\times}$';
            case 'square'
                markers{i} = '\put(0,-0.5){\large{$\blacksquare$}}';
                markers_edge{i} = '\put(-0,-0.5){\large{$\square$}}';
            case 's'
                markers{i} = '\put(0,-0.5){\large{$\blacksquare$}}';
                markers_edge{i} = '\put(-0,-0.5){\large{$\square$}}';
            case 'diamond'
                markers{i} = '\large{$\blacklozenge$}';
                markers_edge{i} = '\large{$\lozenge$}';
            case 'd'
                markers{i} = '\large{$\blacklozenge$}';
                markers_edge{i} = '\large{$\lozenge$}';
            case '^'
                markers{i} = '\large{$\blacktriangle$}';
                markers_edge{i} = '\large{$\vartriangle$}';
            case 'v'
                markers{i} = '\large{$\blacktriangledown$}';
                markers_edge{i} = '\large{$\triangledown$}';
            case '>'
                markers{i} = '\put(0,-0.5){\large{$\blacktriangleright$}}';
                markers_edge{i} = '\put(0,-0.5){\large{$\vartriangleright$}}';
            case '<'
                markers{i} = '\put(0,-0.5){\large{$\blacktriangleleft$}}'; 
                markers_edge{i} = '\put(0,-0.5){\large{$\vartriangleleft$}}';
            case 'pentagram'
                markers{i} = '\large{$\bigstar$}';     
            case 'p'
                markers{i} = '\large{$\bigstar$}';
            case 'hexagram'
                warning('Marker value ''hexagram'' not supported')
                markers{i} = '';
            case 'h'
                warning('Marker value ''h''not supported')
                markers{i} = '';
            otherwise
                warning('Marker value not recognized')
        end
        MarkerFaceColor{i} = get(ud.handles(i), 'MarkerFaceColor');
        if ~isnumeric(MarkerFaceColor{i})
            switch MarkerFaceColor{i}
                case 'auto'
                    MarkerFaceColor{i} = linecolors(i,:);
                case 'none'
                    MarkerFaceColor{i} = [0 0 0];
                    markers{i} = '';
            end
        end
        MarkerEdgeColor{i} = get(ud.handles(i), 'MarkerEdgeColor');
        if ~isnumeric(MarkerEdgeColor{i})
            switch MarkerEdgeColor{i}
                case 'auto'
                    MarkerEdgeColor{i} = linecolors(i,:);
                case 'none'
                    MarkerEdgeColor{i} = [0 0 0];
                    noedge(i) = 1;
                    markers_edge{i} = '';
            end
        end
        if ~isempty(markers{i}) && isempty(markers_edge{i}) && ~noedge(i)
            MarkerFaceColor{i} = MarkerEdgeColor{i};
        end
    else
        warning('No proper style found')
    end
    end
    
    % construct a latex legend string
    finallegendstring = [];
    for i = 1:length(legendstring)
        if ~patch(i)
            finallegendstring_t = [sprintf('\\textcolor[rgb]{%s}{%s} \\put(-19,1){\\textcolor[rgb]{%s}{%s}}  \\put(-19,1){\\textcolor[rgb]{%s}{%s}} & \\hspace{0.5em}%s', regexprep(num2str(linecolors(i,:)),'\s*\s',','), linestyles{i}, regexprep(num2str(MarkerFaceColor{i}),'\s*\s',','), markers{i}, regexprep(num2str(MarkerEdgeColor{i}),'\s*\s',','), markers_edge{i}, legendstring{i})];
            % If the line is a stem, then shift the marker to the end
            if is_stem && ~isempty(markers_edge{i})
                finallegendstring_t = regexprep(finallegendstring_t, 'put\(-19,1\)', 'put\(-5,1\)');
                finallegendstring_t = regexprep(finallegendstring_t, '& \\hspace{0\.5em}', '& \\hspace{1\.5em}');
            end
            finallegendstring = [finallegendstring finallegendstring_t];
        else
            finallegendstring = [finallegendstring sprintf('\\begin{minipage}{3em}\\setlength{\\fboxsep}{0pt}\\setlength{\\fboxrule}{0.5pt}\\fcolorbox[rgb]{%s}{%s}{\\hspace{3em}\\rule{0pt}{1.3ex}}\\end{minipage} & \\hspace{0.5em}%s', regexprep(num2str(patchEdge(i,:)),'\s*\s',','), regexprep(num2str(patchColor(i,:)),'\s*\s',','), legendstring{i})];
        end
        if i ~= length(legendstring)
            finallegendstring = [finallegendstring '\\'];
        end
    end


    
    % calculate the position of the text string
    legendx = min(get(plots(1),'XLim'))+legend_percentual_x*(max(get(plots(1),'XLim')) - min(get(plots(1),'XLim')));
    legendy = min(get(plots(1),'YLim'))+legend_percentual_y*(max(get(plots(1),'YLim')) - min(get(plots(1),'YLim')));

    % print the constructed legend string
    figure(figure_handle)
    axes(plots)
    labels = findobj(figure_handle,'Type','text', 'Tag', 'legend2latex_text', 'Parent', plots(1));
    delete(labels);
    legend_lbl = text(legendx, legendy, sprintf('\\setlength{\\fboxsep}{5pt}\\setlength{\\fboxrule}{0.5pt}\\fcolorbox{black}{white}{\\begin{tabular}{cl} %s \\end{tabular}}', finallegendstring));
    eval(sprintf('set(legend_lbl, %s);', object_opts));
    set(legend_lbl, 'Tag', 'legend2latex_text')

    % delete the legend
%     set(legh,'Visible','Off')
delete(legh)
end %for legh
% Save the Children, so we can restore the original order of axes.
kids = get(figure_handle,'Children');
set(figure_handle,'Children', sort(kids,'descend'));