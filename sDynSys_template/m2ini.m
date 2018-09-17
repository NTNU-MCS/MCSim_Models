function [status, error] = m2ini(System, FileName)
%
% [status, error] = m2ini(System, FileName)
%
% Configuration tool that generates a configuration file with
% all necessary system parameters for the sDynSys S-function model.
%
% Input data:
%
%    System - A Matlab data structure that contains all configuration data for the simulation model, in particular the fields:
%
%       System.SampleTime       Main simulation sample time
%       System.NumContStates    Number of continuous states. Set to 0 for none.
%       System.NumDiscStates    Number of discrete states. Set to 0 for none.
%       System.NumInputSignals  Number of input signals u. Set to 0 for none.
%       System.NumOutputSignals Number of output signals y. Set to 0 for none.
%       System.NumParameters    Number of parameters. Set to 0 for none.
%
%       System.Param{k}.tag  	String with tag of parameter k,     e.g. System.Param{1}.tag = 'A';
%                      .size	Length (int) of parameter vector,   e.g. System.Param{1}.size = 4;
%                      .values  Row vector with space-separated elements (float) of parameter,  
%                                                                   e.g. System.Param{1}.values = [0 1 -2 -1];
%
%                       
%    FileName     - (Required) A string '<filename>', without extension, giving the generated 
%                   configuration file name. 
%
%
% Output data:
%
%    status       - A flag indicating the success of generating the power
%                   configuration data for the S-function implementing a
%                   graph description of the power network. 1
%    error        - text string with an error message. Empty if no error occured.
%
% Generated file:
%    <FileName>.ini - An output file containing all necessary 
%                   configuration data in a well-defined, unambiguous, datastructure 
%                   used by the S-function simuløation model. 
%
%    Copyright: 	Roger Skjetne, NTNU
%    Author:        Roger Skjetne
%    Date created:  2011.04.13  Roger Skjetne.
%    Revised:      	
%
%

%% Initialization
error       = [];   % any possible error message
prjname     = [];   % Project name
bytes       = 0;    % number of bytes written to file
m           = 0;    % counter in datastruct
n           = 0;    % counter in datastruct
p           = 0;
status      = 1;    % Status of execution

% Input argument check
if nargin < 2
    [status, error] = TERMINATE('Not enough input arguments.', 0, fid, filename, bytes);
    return;
elseif nargin > 2
    [status, error] = TERMINATE('Too many input arguments.', 0, fid, filename, bytes);
    return;
end


% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Opening file
% -------------------------------------
filename    = [FileName, '.ini'];
[fid,error] = fopen(filename,'w');

if fid == -1
    error = sprintf('Cannot open specified configuration file: %s. Exiting ...',filename);
    status = -1;
    return;
end

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% Header
% -------------------------------------
if isfield(System,'Project')
    prjname = System.Project;
else
    prjname = '';
end
if isfield(System,'Model')
    model = System.Model;
else
    model = '';
end
count = fprintf(fid,'Simulation model configuration file: %s\n\n',filename); bytes = bytes + count;
count = fprintf(fid,'Project: \t\t%s\n',prjname); bytes = bytes + count;
count = fprintf(fid,'Model:   \t\t%s\n',model); bytes = bytes + count;
count = fprintf(fid,'Date created: \t%s\n',date); bytes = bytes + count;
hh = clock;
count = fprintf(fid,'Time: \t\t\t%s:%s:%s\n\n',num2str(hh(4)),num2str(hh(5)),num2str(hh(6))); bytes = bytes + count;


% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%% System simulation parameters
% -------------------------------------
count = fprintf(fid,'\n'); bytes = bytes + count;
count = fprintf(fid,'System.SampleTime        = %d; \n',System.SampleTime); bytes = bytes + count;
count = fprintf(fid,'System.NumContStates     = %d; \n',System.NumContStates); bytes = bytes + count;
count = fprintf(fid,'System.NumDiscStates     = %d; \n',System.NumDiscStates); bytes = bytes + count;
count = fprintf(fid,'System.NumInputSignals   = %d; \n',System.NumInputSignals); bytes = bytes + count;
count = fprintf(fid,'System.NumOutputSignals  = %d; \n',System.NumOutputSignals); bytes = bytes + count;
count = fprintf(fid,'System.NumParameters     = %d; \n\n',System.NumParameters); bytes = bytes + count;



%% Model parameters
m = 0;
if isfield(System,'Param')
    for i = 1:length(System.Param)
        if ~isempty(System.Param{i})
            flag = isfield(System.Param{i},'values');
            if(flag) 
                ParamSize = size(System.Param{i}.values);
                if(isfield(System.Param{i},'tag'))
                    count = fprintf(fid,'System.Param{%d}.tag   \t = ''%s''; \n', m, System.Param{i}.tag); bytes = bytes + count;
                else
                    count = fprintf(fid,'System.Param{%d}.tag   \t = NoTag; \n', m); bytes = bytes + count;
                end
                count = fprintf(fid,'System.Param{%d}.size  \t = %d; \n', m, prod(ParamSize)); bytes = bytes + count;
                count = fprintf(fid,'System.Param{%d}.values\t = [', m); bytes = bytes + count;
                for j=1:ParamSize(1)
                    for k=1:ParamSize(2)
                        count = fprintf(fid,'%.3e ', System.Param{i}.values(j,k)); bytes = bytes + count;
                    end
                end
                count = fprintf(fid,']; \n\n'); bytes = bytes + count;
%                 count = fprintf(fid,'%.3e]; \n\n', System.Param{i}.values(System.Param{i}.size)); bytes = bytes + count;
                m = m + 1;
            end
        end
    end
else 
    [status, error] = TERMINATE(sprintf('Configuration data not contained in System. Cannot complete operation. '), 0, fid, filename, bytes);
    return;
end


%% Terminating function
[status, error] = TERMINATE('', 1, fid, filename, bytes);






%% Subfunction: TERMINATE
%  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
function [status, error] = TERMINATE(error_msg,status_flag,file_id, file_name,num_bytes)

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% Closing file
% -------------------------------------
flag = fclose(file_id); 

if flag == -1
    fprintf(1,'Cannot close file: %s\n', file_name);
    status = 0;
    error  = strcat(error_msg, '\n', sprintf('Cannot close file: %s\n', file_name));
elseif status_flag == 0
    fprintf(1,'Error: %s\n File closed.\n',error_msg);
    status = 0;
    error  = error_msg;
elseif status_flag == 1
    fprintf(1,'%d bytes written to file: %s\n',num_bytes,file_name);
    status = status_flag;
    error  = error_msg;
else
    fprintf(1,'An unknown error has occurred.');
    status = 0;
    error  = [error_msg, '\nAn unknown error has occurred.\n'];
end


