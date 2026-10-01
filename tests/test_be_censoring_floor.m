% Floor mode of be_censoring: FAIL rows and NaN data errors are never washed; backward-error
% mode wins when the header carries a positive be_tol.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'plotting'));
tmp = [tempname, '.csv'];
fid = fopen(tmp, 'w'); fprintf(fid, '# synthetic noise_level=1e-11\nalgorithm,m,n\nx,75824,8304\n'); fclose(fid);
fl = 1e-11 * sqrt(1 - 8304/75824);
T = table({'at_floor';'above';'above_failed';'no_value'}, repmat(75824,4,1), repmat(8304,4,1), ...
          [fl; 3*fl; 3*fl; NaN], [0; 0; 1; 0], ...
          'VariableNames', {'algorithm','m','n','ls_data_error','qr_status'});
C = be_censoring(T, (1:4).', T.qr_status ~= 0, tmp);
assert(strcmp(C.mode, 'floor'));
assert(abs(C.floor - fl) < 1e-20);
assert(isequal(C.censored,  [false; true; false; false]));
assert(isequal(C.converged, [true; false; false; false]));
assert(~any(C.no_oracle));
% precedence: a positive be_tol in the header selects backward-error mode
fid = fopen(tmp, 'w'); fprintf(fid, '# be_tol=3.05e-15 noise_level=1e-11\nalgorithm\nx\n'); fclose(fid);
T2 = table({'a';'b'}, [1;1], [1;1], [1;1], [0;0], [1e-16; 1e-14], {'be';'floor'}, [1;1], ...
           'VariableNames', {'algorithm','m','n','ls_data_error','qr_status','be_final','stop_reason','t_be_us'});
C2 = be_censoring(T2, (1:2).', false(2,1), tmp);
assert(strcmp(C2.mode, 'be'));
assert(isequal(C2.censored, [false; true]));
delete(tmp);
disp('test_be_censoring_floor: OK');
