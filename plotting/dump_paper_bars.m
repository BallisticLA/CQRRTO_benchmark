function dump_paper_bars(repo, out_json)
% Run the paper plotters for era 1002_p5 on invisible figures and dump every
% plotted bar value, tick label, axis label and text annotation to JSON.
% Cell discovery and CSV selection copy run_all.m (newest CSV, merge_gram_arms).
set(groot, 'DefaultFigureVisible', 'off');
addpath(repo); addpath(fullfile(repo, 'utils')); addpath(fullfile(repo, 'plotting'));
out = struct('cells', {{}});

toep_dir = fullfile(repo, 'results', 'toeplitz_ls_1002_p5_pcg_ne');
commit_t = assert_era_commit(toep_dir, '[1002 p5 in-loop-BE]');
for s = {'fixedm', 'small', 'middle', 'large'}
    size_dir = fullfile(toep_dir, s{1});
    merged = merge_gram_arms(size_dir, '*_toeplitz_ls_results.csv');
    if ~isempty(merged)
        size_dir = fullfile(size_dir, 'merged'); res_csv = merged;
    else
        f = dir(fullfile(size_dir, '*_toeplitz_ls_results.csv')); f = f([f.bytes] > 0);
        [~, ord] = sort({f.name}); res_csv = f(ord(end)).name;
    end
    close all;
    plot_toeplitz_results(size_dir, res_csv, s{1}, [], 'best');
    c = harvest(gcf); c.bench = 'toeplitz'; c.cell = s{1}; c.csv = fullfile(size_dir, res_csv); c.commit = commit_t;
    out.cells{end+1} = c;
end

camp_dir = fullfile(repo, 'results', 'irlsq_reg_1002_p5');
commit_f = assert_era_commit(camp_dir, '[1002 p5 in-loop-BE]');
for s = {'native_ill_dd', 'small_kc1e10_dd', 'large_kc1e10_dd'}
    cell_dir = fullfile(camp_dir, s{1});
    merged = merge_gram_arms(cell_dir, '*_irlsq_reg_results.csv');
    if ~isempty(merged)
        cell_dir = fullfile(cell_dir, 'merged'); res_csv = merged;
    else
        f = dir(fullfile(cell_dir, '*_irlsq_reg_results.csv')); f = f([f.bytes] > 0);
        [~, ord] = sort({f.name}); res_csv = f(ord(end)).name;
    end
    bd_csv = strrep(res_csv, '_results.csv', '_breakdown.csv');
    close all;
    plot_irlsq_results(cell_dir, res_csv, bd_csv, s{1}, [], [], 'best');
    c = harvest(gcf); c.bench = 'fem'; c.cell = s{1}; c.csv = fullfile(cell_dir, res_csv); c.commit = commit_f;
    out.cells{end+1} = c;
end

fid = fopen(out_json, 'w'); fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true)); fclose(fid);
fprintf('wrote %s (%d cells)\n', out_json, numel(out.cells));
end

function c = harvest(fig)
axs = findobj(fig, 'Type', 'axes');
tiles = arrayfun(@(a) a.Layout.Tile, axs);
[~, o] = sort(tiles); axs = axs(o);
c = struct('panels', {{}});
for a = axs(:)'
    p = struct();
    p.tile = a.Layout.Tile;
    p.title = char(string(a.Title.String));
    p.ylabel = char(string(a.YLabel.String));
    p.yscale = a.YScale; p.ylim = a.YLim;
    p.xticklabels = cellstr(a.XTickLabel);
    bars = findobj(a, 'Type', 'bar');
    bars = flipud(bars);   % creation order
    p.bars = {};
    for b = bars(:)'
        q = struct('x', b.XData, 'y', b.YData, 'facecolor', b.FaceColor, 'facealpha', b.FaceAlpha, ...
                   'linestyle', b.LineStyle, 'edgecolor', b.EdgeColor, 'name', b.DisplayName);
        p.bars{end+1} = q;
    end
    txt = findobj(a, 'Type', 'text');
    p.texts = {};
    for t = txt(:)'
        if isempty(t.String), continue; end
        p.texts{end+1} = struct('s', char(string(t.String)), 'pos', t.Position);
    end
    c.panels{end+1} = p;
end
end
