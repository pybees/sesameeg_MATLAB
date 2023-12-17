function varargout = process_posterior( varargin )
% PROCESS_POSTERIOR: 

% @=============================================================================
% This function is part of the SESAME software:
%
% 
% Copyright (c) 
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF XXX AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "xxx license" at command prompt.
% =============================================================================@
%
% Authors: 

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Plot sources';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 1001;
    sProcess.Description = 'Plot sources';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some options
    sProcess.isSourceAbsolute = -1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    ProtocolInfo = bst_get('ProtocolInfo');
    sStudy = bst_get('Study', ProtocolInfo.iStudy);

    ResultsMat = load(fullfile(ProtocolInfo.STUDIES, sInput(1).FileName));

    vertices = ResultsMat.GridLoc;
    posterior = ResultsMat.posterior;
    
    % === NOT SHARED ===
    % Get the output sStudy (pick the one from the first file)
    
    iStudy = sInput(1).iStudy;
    % Create a default output filename

    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInput(1).FileName), 'dipoles')

    sesame_dipoles = db_template('dipolemat');
    sesame_dipoles.Comment = 'Sesame est dipoles ';
    sesame_dipoles.Time = 0;
    sesame_dipoles.DipoleNames = 'estimated dipoles';
    sesame_dipoles.DataFile = sInput(1).FileName;

    n_est_dip = length(posterior.estimated_dipoles);
    nt = length(ResultsMat.cfg.time_window);
    for nd = 1:n_est_dip
      sesame_dipoles.Dipole(nd).Index = 1;
      sesame_dipoles.Dipole(nd).Time = 0;
      sesame_dipoles.Dipole(nd).Origin = [0 0 0];
      sesame_dipoles.Dipole(nd).Loc = vertices(posterior.estimated_dipoles(nd), :)';
      sesame_dipoles.Dipole(nd).Amplitude = posterior.QV_estimated(3*nd-2:3*nd, end);
      sesame_dipoles.Dipole(nd).Goodness = posterior.gof;  
    end
    
    % Add to sStudy
    iDipole = length(sStudy.Dipoles) + 1;
    sStudy.Dipoles(iDipole).FileName = file_short(OutputFiles{1});
    sStudy.Dipoles(iDipole).Comment  =  'estimated dipoles ';
    sStudy.Dipoles(iDipole).DataFile = sesame_dipoles.DataFile;
    
    bst_set('Study', iStudy, sStudy);
    % Save sStudy
    bst_save(OutputFiles{1}, sesame_dipoles)
    db_add_data(iStudy, OutputFiles{1}, sesame_dipoles);
    % Update tree
    panel_protocols('UpdateNode', 'sStudy', iStudy);
end




