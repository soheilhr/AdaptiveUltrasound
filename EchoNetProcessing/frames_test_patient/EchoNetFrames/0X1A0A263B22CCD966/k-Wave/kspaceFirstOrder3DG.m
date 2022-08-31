function sensor_data = kspaceFirstOrder3DG(varargin)
%KSPACEFIRSTORDER3DG 3D time-domain simulation of wave propagation on a GPU using C++ CUDA code.
%
% DESCRIPTION:
%     kspaceFirstOrder3DG provides a blind interface to the C++/CUDA
%     version of kspaceFirstOrder3D (called kspaceFirstOrder-CUDA) in the
%     same way as kspaceFirstOrder3DC. Note, the C++ code does not support
%     all input options, and all display options are ignored (only command
%     line outputs are given). See the k-Wave user manual for more
%     information. 
%
%     The function works by appending the optional input 'SaveToDisk' to
%     the user inputs and then calling kspaceFirstOrder3D to save the input
%     files to disk. The contents of sensor.record (if set) are parsed as
%     input flags, and the C++ code is run using the system command. The
%     output files are then automatically loaded from disk and returned in
%     the same fashion as kspaceFirstOrder3D. The input and output files
%     are saved to the temporary directory native to the operating system,
%     and are deleted after the function runs.
%
%     This function requires the C++ binary/executable of 
%     kspaceFirstOrder-CUDA to be downloaded from
%     http://www.k-wave.org/download.php and placed in the "binaries"
%     directory of the k-Wave toolbox. Alternatively, the name and location
%     of the binary can be specified using the optional input parameters
%     'BinaryName' and 'BinariesPath'.
%
% USAGE:
%     see kspaceFirstOrder3D, kspaceFirstOrder3DC
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 30th October 2013
%     last update     - 6rd March 2019
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2019 Bradley Treeby
%
% See also kspaceFirstOrder3D, kspaceFirstOrder3DC

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% This function is essentially a wrapper and directly uses the capabilities
% of kspaceFirstOrder3DC by replacing the binary name with the name of the
% GPU binary. 

% Check for the binary name input. If not defined, set the default name of
% the GPU binary
if ~any(strcmp('BinaryName', varargin))
    if isunix
        varargin = [varargin {'BinaryName', 'kspaceFirstOrder-CUDA'}];
    else
        varargin = [varargin {'BinaryName', 'kspaceFirstOrder-CUDA.exe'}];
    end
end

% pass inputs to kspaceFirstOrder3DC
sensor_data = kspaceFirstOrder3DC(varargin{:});