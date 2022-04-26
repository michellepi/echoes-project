function [b,a]=make_digital_filter(f_cut,fs,parameters);
% calculate the parameters a and b for a digital filter
% f_cut is the cut-off frequency, paramters.lphp='low','high','stop' (in the latter, f_cut must hold both cut-off frequencies)
% parameters is an optionl input (see default values below)
% and parameters.type='butter','cheby1', 'cheby2', or 'ellip'
% for 'butter', parameters.N gives the order
% for 'cheby1', parameters.N gives the order, and parameters.R gives the ripple (typically R=0.5)
% for 'cheby2', parameters.N gives the order, and parameters.R gives the ripple (typically R=20)
% for 'ellip', parameters.N gives the order, parameters.Rp and parameters.Rs give the ripple in the pass and stop-bands respectively (typial 0.5 and 20)

if nargin==2
    parameters.N=4;
end
if ~isfield(parameters,'N') % set default filter order to 4
    parameters.N=4;
end
if ~isfield(parameters,'type'); % set default filter type to Butterworth
    parameters.type='butter';
end
if ~isfield(parameters,'lphp'); % set default filter type to low-pass
    parameters.lphp='low';
end
if~isfield(parameters,'R');
    parameters.R=0.5;
end
if~isfield(parameters,'Rp');
    parameters.R=0.5;
end
if~isfield(parameters,'Rs');
    parameters.R=20;
end
switch lower(parameters.type);
case('butter');
    [b,a]=butter(parameters.N,2*f_cut/fs,parameters.lphp);
case('cheby1');
    [b,a]=cheby1(parameters.N,parameters.R,2*f_cut/fs,parameters.lphp);
case('cheby2');
    [b,a]=cheby2(parameters.N,parameters.R,2*f_cut/fs,parameters.lphp);
case('ellip');
    [b,a]=ellip(parameters.N,parameters.Rp,parameters.Rs,2*f_cut/fs,parameters.lphp);
otherwise
    disp('ERROR IN MAKE_DIGITAL_FILTER: No known filter type');
end

