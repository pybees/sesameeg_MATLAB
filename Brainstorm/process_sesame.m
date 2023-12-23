function varargout = process_sesame( varargin )
% PROCESS_SESAME: Call SESAME algorithm
% 
% inverse_SESAME samples the posterior distribution of a
% multi-dipole Bayesian model and provides an estimate of the number of
% dipoles and of the dipole locations and time courses; in addition, it
% provides uncertainty quantification in terms of a posterior probability
% map. Dipoles are assumed to have fixed location during the analyzed time
% window.
% 
% REFERENCES: 
% 
% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as publishebead by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: 

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description of the process
    sProcess.Comment     = 'SESAME';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';  % sub-menu in which the process appears
    sProcess.Index       = 1000;
    sProcess.Description = 'https://pybees.github.io/sesameeg_MATLAB';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;  % appears in Process1 tab
    sProcess.nMinFiles   = 1;

    % Definition of the options
    % Options: Time window
    sProcess.options.label1.Comment = '<B>Input options</B>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.active_time.Comment = 'Time window: ';
    sProcess.options.active_time.Type    = 'timewindow';
    sProcess.options.active_time.Value   = [];
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    sProcess.options.sensortypes.InputTypes = {'data'};

    % Options: SESAME parameters
    sProcess.options.label2.Comment = '<BR><B>SESAME options</B>:';
    sProcess.options.label2.Type    = 'label';
    % Options: number of samples
    sProcess.options.particles.Comment = 'Number of particles: ';
    sProcess.options.particles.Type    = 'value';
%     sProcess.options.particles.Class   = 'Advaced_options';
    sProcess.options.particles.Value   = {100, '', 0};
    % Options: noise std
    sProcess.options.noise_std.Comment = 'noise standard deviation(optional):';
    sProcess.options.noise_std.Type    = 'value';
    sProcess.options.noise_std.Value   = {[]};    
    % Options: dip mom std
    sProcess.options.dip_mom_std.Comment = 'dipole moment standard deviation(optional):';
    sProcess.options.dip_mom_std.Type    = 'value';
    sProcess.options.dip_mom_std.Value   = {[]};
%     sProcess.options.dip_mom_std.Class   = 'Advaced_options';
%     sProcess.options.label3.Comment = '<BR><B>Advanced options</B>:';
%     sProcess.options.label3.Type    = 'label';
%     sProcess.options.label3.Controller = 'Advanced_options'; %  = 'Advanced options';
    % Options: hyper q
    sProcess.options.hyper_q.Comment = 'Use hyperprior';
    sProcess.options.hyper_q.Type    = 'checkbox';
%     sProcess.options.hyper_q.Class   = 'Advaced_options';
    sProcess.options.hyper_q.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned list of files
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    ActiveTime  = sProcess.options.active_time.Value{1}
    SensorTypes = sProcess.options.sensortypes.Value
    cfg.noise_std     = sProcess.options.noise_std.Value{1}
    cfg.dip_mom_std     = sProcess.options.dip_mom_std.Value{1}
    cfg.bool_hyper_q     = sProcess.options.hyper_q.Value
    cfg.n_samples   = sProcess.options.particles.Value{1}
    
    % ===== LOAD CHANNEL FILE =====
    % Load channel file
    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    % Find the MEG channels
%     iMEG = good_channel(ChannelMat.Channel, [], 'MEG')
%     iEEG = good_channel(ChannelMat.Channel, [], 'EEG');
    iChannels = channel_find(ChannelMat.Channel, SensorTypes);
    
    % ===== LOAD HEAD MODEL =====
    % Get channel study
    [sChannelStudy, ~] = bst_get('ChannelFile', sInputs(1).ChannelFile);
    % Load the default head model
    HeadModelFile = sChannelStudy.HeadModel(sChannelStudy.iHeadModel).FileName;
    sHeadModel = load(file_fullpath(HeadModelFile));
    % Get number of sources
    nSources = length(sHeadModel.GridLoc);
    
    % ===== LOAD THE DATA =====
    DataMat = in_bst(sInputs(1).FileName, [], 0);
    sInputs(1).FileName
    nChannels = size(DataMat.F,1);
    nTime     = size(DataMat.F,2);
    Time      = DataMat.Time;

    if isempty(ActiveTime)
      iActiveTime = 1:nTime;
    else
      iActiveTime = panel_time('GetTimeIndices', Time, ActiveTime);
    end
    cfg.t_start = iActiveTime(1);
    cfg.t_stop = iActiveTime(end);
    cfg.time_window = Time(cfg.t_start:cfg.t_stop)

    % Remove bad channels
    iBadChan = find(DataMat.ChannelFlag == -1);
    iChannels = setdiff(iChannels, iBadChan);
    % Error: All channels tagged as bad
    if isempty(iChannels)
        bst_report('Error', sProcess, sInputs, 'All the selected channels are tagged as bad.');
        return;
    end
    
    % ===== MODEL SCALING =====
    % Model scaling, definition of the scaling parameter theta^*
    LF = sHeadModel.Gain(iChannels, :);
    B = DataMat.F(iChannels, :);
    size(B)
    size(LF)
    vertices = sHeadModel.GridLoc;

    % ===== PROCESS =====
    % Run the IAS algorithm
    disp('Running SESAME ...');
    posterior = Compute(B, LF, vertices, cfg);
    
    n_est_dip = length(posterior.estimated_dipoles)
    ImageGridAmp = sum(posterior.pmap(:,posterior.final_it, n_est_dip),3);

    % ===== SAVE THE RESULTS =====
    % Create a new data file structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = ImageGridAmp;
%     ResultsMat.nComponents   = 3;   % 1 (constrained) or 3 (uncostrained)
    ResultsMat.Comment       = ['Sesame_pmap' set_comment(cfg, Time)];
    ResultsMat.Function      = 'Sesame';
    ResultsMat.Time          = ActiveTime;  % Leave it empty if using ImagingKernel
    ResultsMat.DataFile      = sInputs(1).FileName;
    ResultsMat.HeadModelFile = HeadModelFile;
    ResultsMat.HeadModelType = sHeadModel.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = iChannels;
    ResultsMat.SurfaceFile   = sHeadModel.SurfaceFile;
    ResultsMat.GridLoc       = vertices;
    ResultsMat.cfg           = cfg;
    ResultsMat.posterior     = posterior;

    % === NOT SHARED ===
    % Get the output study (pick the one from the first file)
    [sStudy, iStudy] = bst_get('Study', sInputs.iStudy);
    % Create a default output filename 
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'results')
    % Save on disk
    bst_save(OutputFiles{1}, ResultsMat)
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, ResultsMat);
    % Update the database explorer display
    panel_protocols('UpdateNode', 'Study', iStudy);

end


%% ===== EXTERNAL CALL =====
% USAGE: posterior = process_sesame('Compute', B, LF, vertices, cfg)
function posterior = Compute(B, LF, vertices, cfg)
  posterior = inverse_SESAME(B, LF, vertices, cfg);
end


function str_ris = set_comment(cfg, Time)

  if cfg.bool_hyper_q == 0
    str_hq = '_no_hq';
  else
    str_hq = '';
  end

  if cfg.t_start == cfg.t_stop
    str_t = sprintf('t_[%0.*f]ms', 1, Time(cfg.t_start));
  else
    str_t = [sprintf('tw_[%0.*f', 1, Time(cfg.t_start)*1000), '_', sprintf('%0.*f]ms', 1, Time(cfg.t_stop)*1000)];
  end
  
  str_ris = [str_hq, '_', str_t]

end

