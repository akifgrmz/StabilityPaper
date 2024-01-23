function msfuntmpl_basic(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;
block.NumDialogPrms  = 4;

tstep = block.DialogPrm(3).Data;
block.SampleTimes = [-1 0];
FrameLength=block.DialogPrm(2).Data;

% Setup port properties to be inherited or dynamic
% block.SetPreCompInpPortInfoToDynamic;
% block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = FrameLength;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;
% 
% block.InputPort(2).Dimensions        = 1;
% block.InputPort(2).DatatypeID  = 0;  % double
% block.InputPort(2).Complexity  = 'Real';
% block.InputPort(2).DirectFeedthrough = true;

% Override output port properties
block.OutputPort(1).Dimensions  = FrameLength;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';
% 
% block.OutputPort(2).Dimensions       = 1;
% block.OutputPort(2).DatatypeID  = 0; % double
% block.OutputPort(2).Complexity  = 'Real';
% Register parameters

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)

block.NumDworks = 3;
FrameLength = block.DialogPrm(2).Data;
GSOrder = block.DialogPrm(1).Data;


  %dwork 1 
%       label=sprintf('x%d',OrderI);
    block.Dwork(1).Name            = 'x1';
    block.Dwork(1).Dimensions      = FrameLength*(GSOrder+1);
    block.Dwork(1).DatatypeID      = 0;      % double
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

  %dwork 2 
  block.Dwork(2).Name            = 'x2';
  block.Dwork(2).Dimensions      = FrameLength;
  block.Dwork(2).DatatypeID      = 0;      % double
  block.Dwork(2).Complexity      = 'Real'; % real
  block.Dwork(2).UsedAsDiscState = true;
%   
  %dwork 3 
  block.Dwork(3).Name            = 'x3';
  block.Dwork(3).Dimensions      = 1;
  block.Dwork(3).DatatypeID      = 0;      % double
  block.Dwork(3).Complexity      = 'Real'; % real
  block.Dwork(3).UsedAsDiscState = true;
%   
%   %dwork 4 
%   block.Dwork(4).Name            = 'x4';
%   block.Dwork(4).Dimensions      = 1;
%   block.Dwork(4).DatatypeID      = 0;      % double
%   block.Dwork(4).Complexity      = 'Real'; % real
%   block.Dwork(4).UsedAsDiscState = true;
%   
%   %dwork 5 
%   block.Dwork(5).Name            = 'x5';
%   block.Dwork(5).Dimensions      = 1;
%   block.Dwork(5).DatatypeID      = 0;      % double
%   block.Dwork(5).Complexity      = 'Real'; % real
%   block.Dwork(5).UsedAsDiscState = true;
% 
%    %dwork 6
%   block.Dwork(6).Name            = 'x6';
%   block.Dwork(6).Dimensions      = 1;
%   block.Dwork(6).DatatypeID      = 0;      % double
%   block.Dwork(6).Complexity      = 'Real'; % real
%   block.Dwork(6).UsedAsDiscState = true;

%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)

%  block.Dwork(1).Data = zeros(100,1);
% block.Dwork(2).Data = zeros(100,1);
% %end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)
%collect data

FrameLength = block.DialogPrm(2).Data;
M=block.DialogPrm(1).Data;
in=block.InputPort(1).Data;
store=zeros((M+1)*FrameLength,1);

store(1:end-FrameLength)=block.Dwork(1).Data(FrameLength+1:end);
store(end-FrameLength+1:end)=in;
block.Dwork(1).Data=store;

EMGFrame=reshape(store,M+1,FrameLength);
[gsFrame,~]=gsfm(EMGFrame,M);

block.OutputPort(1).Data=gsFrame;


% end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)




%block.Dwork(1).Data = block.InputPort(1).Data;

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

function [y_out, m_out] = gsfm(x, M)
% general Gram-Schmidt Prediction Error Filter
% adapted from Nathan Makowsky's 'SimpleGSFilt.m'
% x is M+1 frames of row vectors of length n_f
% frames are in reverse temporal order - current frame is index 1

[~, n_per_frame] = size(x);
w = zeros(M);
eps = zeros(M+1, M+1, n_per_frame);
for i = 1:M+1
    eps(i, 1, :) = x(i,:);
end
for m = 1:M % step m = m_step-1 through 0:M-1
    for i = 1:1:(M-m+1) % step i = i_step-1 through 0:M-m-1
%         eps(i_step, m_step, :) = x(m_step,i_step,:);
        temp1(1, :) = eps(i, m, :);
        temp2(:, 1) = eps(M-m+2, m, :);
        w(m, i) = temp1*temp2/(mag(eps(M-m+2, m, :)))^2;
        eps(i, m+1, :) = eps(i, m, :) - w(m, i)*eps(M-m+2, m, :);
    end
end
y_out(1, :) = eps(1, M+1, :);
m_out(1, :) = x(1, :) - y_out(1, :);
 % of function gsfm
function y = mag(X)
y = sqrt(sum((X.^2)));
