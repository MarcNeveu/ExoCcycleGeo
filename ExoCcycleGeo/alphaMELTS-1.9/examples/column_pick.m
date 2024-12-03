function [data, names] = column_pick(column_list_file, varargin)
%COLUMN_PICK run column_pick.command and process results.
%
%   [DATA, NAMES] = COLUMN_PICK(COLUMN_LIST_FILE, ...)
%
%   Calls column_pick.command for the given COLUMN_LIST_FILE using Matlab's 
%   perl function. Additional COLUMN_LIST_FILES can follow the first and
%   will behave as if concatenated.
%
%   The output is a matrix DATA, padded with NaN where appropriate, and
%   (optionally) a cell NAMES with the column_pick.command generated column
%   names. You can assign the columns of the data matrix to variabls with
%   these names (or your own choice) with something like this:
%
%   for k = 1:length(names)
%       eval([names{k}, '= data(:, ', int2str(k), ');']);
%   end
%
%   COLUMN_PICK assumes that the COLUMN_LIST_FILES uses the following lines:
%
%   Delimiter: ','
%   Header: matlab
%
%   It should be obvious if either the delimiter or header is incorrect. If
%   cell2mat fails it could mean there is a misake in the COLUMN_LIST_FILE
%   (e.g. 'variables' were not used in a table with 'missing' rows) or there
%   is a bug in column_pick. Try running perl('column_pick.command', ...)
%   without the '> ' to see if the script generated any warnings or errors.   
%
%   See alphaMELTS documentation and forum for more information. We have 
%   thoroughly tested column_pick.command and believe it is robust but if 
%   you find a problem please report it (psmith@gps.caltech.edu). Thanks!

if (nargin < 1); error('Please supply a column_list_file name!'); end;

data = perl('column_pick.command', column_list_file, varargin{:}, '> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

data = textscan(data, '%s', 'Delimiter', delim);
data = data{:}';

names = textscan(data{1}, '%s', 'Delimiter', ',');
names = names{:}';

data = cellfun(@(x) textscan(x, '%f', 'Delimiter', ','), data(2:end));
data = cell2mat(data)';

end
