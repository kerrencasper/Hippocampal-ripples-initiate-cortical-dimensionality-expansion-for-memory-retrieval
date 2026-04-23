function paths = config_paths()

% === USER-DEFINED ROOT PATH ===
paths.base_path = '/path_to_downloaded_data_from_Zenodo/'; % Note! Not possible to share raw/processed data, only data to reproduce figures

% === PROJECT STRUCTURE ===
paths.data_dir          = fullfile(paths.base_path, 'preprocessing', 'artifact_rejected_data');
paths.save_dir          = fullfile(paths.base_path, 'output_data', 'decoding');
paths.channels_dir      = fullfile(paths.base_path, 'revised_data', 'templates');
paths.anatomy_dir       = fullfile(paths.base_path, 'revised_data', 'additional_analyses', 'visualisation');
paths.AAL_dir           = fullfile(paths.base_path, 'revised_data', 'subfunctions', 'AAL3');
paths.SPM_dir           = fullfile(paths.base_path, 'revised_data', 'subfunctions', 'spm12');

% === subfunctions ===

paths.add_subfunctions = '/path_to_downloaded_data_from_Zenodo/';

paths.help_functions = fullfile(paths.add_subfunctions, 'revised_data', 'help_functions');
paths.MVPA_Light_master = fullfile(paths.add_subfunctions, 'revised_data', 'MMVPA_Light_master');
paths.plotting = fullfile(paths.add_subfunctions, 'revised_data', 'plotting');
paths.subfunctions = fullfile(paths.add_subfunctions, 'revised_data', 'subfunctions');


end