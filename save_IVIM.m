function [] = save_IVIM(volume,filename_for_header,varargin)
%SAVE_IVIM  Save IVIM parameter maps to disk (4D NIfTI or text file).
%
%   SAVE_IVIM(VOLUME, FILENAME_FOR_HEADER, ...) saves the fitted IVIM
%   parameter maps (S0, f, D*, D) using the header information from an
%   existing NIfTI file.
%
%   Required inputs
%   ---------------
%   VOLUME              : 4D numeric array
%       IVIM parameter maps to save, with volumes ordered as:
%         1. S0 (baseline signal)
%         2. f  (perfusion fraction)
%         3. D* (pseudo-diffusion coefficient)
%         4. D  (diffusion coefficient)
%
%   FILENAME_FOR_HEADER : char/string
%       Path to an existing NIfTI file used as a header template
%       (typically the original DWI data). Provides orientation,
%       voxel size, and other metadata for saving.
%
%   Name–Value pairs (optional)
%   ---------------------------
%   'output_filename' : char/string (default: FILENAME_FOR_HEADER)
%       Base filename for the output (extension is added automatically).
%
%   'output_path'     : char/string (default: "./")
%       Directory where the output file will be saved.
%
%   'prefix'          : char/string (default: "IVIM_")
%       Prefix for the saved result filename.
%
%   Output
%   ------
%   None (file written to disk).
%
%   Notes
%   -----
%   - If VOLUME is a vector of length 4, results are saved as a text file
%     with extension `.txt`.
%   - Otherwise, results are saved as a 4D NIfTI file with dimensions
%     matching the header template.
%   - Output NIfTI will always have 4 volumes in the 4th dimension
%     (S0, f, D*, D).
%
%   Example
%   -------
%     save_IVIM(results, 'sub01_dwi.nii', 'prefix', 'fit_', ...
%               'output_path', './results/');
%
p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'volume');
addRequired(p,'filename_for_header'); % file to read template header (like original dwi data) and determine size and orientation to save
addParameter(p,'output_filename',filename_for_header);
addParameter(p,'output_path',"");
addParameter(p,'prefix','IVIM_');

parse(p,volume,filename_for_header,varargin{:});


if contains(filename_for_header,"/")
    newname=strsplit(filename_for_header,"/");
    output_path=strjoin(newname(1:end-1),"/")+"/";
    name=newname(end);
else
    output_path="./";
    name=p.Results.output_filename;
end
if numel(volume)==4
    writematrix(volume, output_path+p.Results.prefix+name+".txt");
else
    header=niftiinfo(filename_for_header);
    header.ImageSize(4)=4;
    header.Datatype=class(volume);
    niftiwrite(volume,output_path+p.Results.prefix+name,header)

end

end