function RS=restingstatenoise(fs,fo,bw,T,vx)
% RESTINGSTATENOISE generates vx voxels of T samples at fs Hz, with alpha
% components characterized by central freq. fo and banwidth bw, in Hz
%
%   RS=restingstatenoise(fs,fo,bw,T,vx)
%
%   fs = sampling freq in Hz (default = 200)
%   fo = central freq in Hz of alpha component (default = 10 Hz)
%   bw = bandwidth in Hz of alpha component (default = 2 Hz)
%   T = number of samples per channel (default = 201, i.e. 1 s)
%   vx = number of voxels (default = 6003)
%   RS = resting state noise in the voxels

% atc 2 deirel, 2017-03-09

%% default inputs
if (nargin<5)||isempty(vx)
    vx=6003;
end
if (nargin<4)||isempty(T)
    T=201;
end
if (nargin<3)||isempty(bw)
    bw=2;
end
if (nargin<2)||isempty(fo)
    fo=10;
end
if (nargin<1)||isempty(fs)
    fs=200;
end
%%
RS=randn(vx,T);
[b,a] = iirpeak(fo/(fs/2),bw/(fs/2));   % Second-order IIR peak or resonator filter to generate alpha
%fvtool(b,a);
tmp=filter(b,a,RS');
RS=tmp';
%% 

