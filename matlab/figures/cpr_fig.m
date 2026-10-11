function S = cpr_fig()
%CPR_FIG  Figure conventions for the paper, as a struct of helpers.
%Template: matlab/scripts/make_corrected_figs.m (Figs. 15-18 of the old paper):
%  12 x 6 cm page, print -dpdf -painters (vector PDF), included in LaTeX with
%  height=0.5\columnwidth; LineWidth 1.5; no markers; line styles '-','--','-.',':';
%  MATLAB default colours with YELLOW third; explicit legend locations (never 'best', which
%  silently drops the legend on invisible figures printed with -painters); grid and box on.
%  Font: Times New Roman, default size. Figures stay open on screen and are also written to PDF and PNG when
%  cpr_settings().figures.print is true. cpr_settings().figures.axes selects the paper's axis limits or automatic ones.
%   S = cpr_fig();  f = S.new();  ...plot...;  S.style(gca);  S.save(f, file)   (12 x 6 cm)
%   f = S.new(24, 6); ... S.save(f, file, 24, 6)                                 (other sizes)
S.col = struct('IS', [0 0.447 0.741], 'RS', [0.850 0.325 0.098], 'IC', [0.929 0.694 0.125], 'RC', [0.494 0.184 0.556]);
S.cols = [S.col.IS; S.col.RS; S.col.IC; S.col.RC; 0.466 0.674 0.188];
S.ls = struct('IS', '-', 'RS', '--', 'IC', '-.', 'RC', ':');
S.lss = {'-', '--', '-.', ':', '-'};
S.lw = 1.5;
S.font = 'Times New Roman';
S.new = @newfig;
S.style = @styleax;
S.save = @savefig;
S.paperAxes = @() strcmp(cpr_settings().figures.axes, 'paper');   % true: the paper's axis limits and ticks; false: automatic
end

function f = newfig(w, h)
if nargin < 1, w = 12; h = 6; end
f = figure('Visible', 'on', 'Renderer', 'painters', 'Units', 'centimeters', 'Position', [2 2 w h], 'Color', 'w');
fn = 'Times New Roman';
set(f, 'DefaultAxesFontName', fn, 'DefaultTextFontName', fn, 'DefaultLegendFontName', fn, ...
    'DefaultAxesLineWidth', 0.5, 'DefaultLineLineWidth', 1.5, ...
    'DefaultAxesFontSize', 11, 'DefaultLegendFontSize', 11);   % 11 pt
end

function styleax(ax)
grid(ax, 'on'); box(ax, 'on'); set(ax, 'FontName', 'Times New Roman', 'Layer', 'top');
end

function savefig(f, file, w, h)
if nargin < 3, w = 12; h = 6; end
set(f, 'PaperUnits', 'centimeters', 'PaperSize', [w h], 'PaperPosition', [0 0 w h]);
[d, nm] = fileparts(file); if isempty(d), d = pwd; end
if cpr_settings().figures.print
    for attempt = 1:3                  % printing occasionally fails when other MATLAB sessions use graphics
        try
            print(f, fullfile(d, [nm '.pdf']), '-dpdf', '-painters');
            print(f, fullfile(d, [nm '.png']), '-dpng', '-r150');
            break
        catch err
            if attempt == 3, rethrow(err); end
            pause(3);
        end
    end
end
drawnow;
end
