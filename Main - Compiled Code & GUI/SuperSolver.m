% This GUI is based on using Simple Optimizer Tabs from the WFAToolbox Team
function varargout = SuperSolver(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SuperSolver_OpeningFcn, ...
    'gui_OutputFcn',  @SuperSolver_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SuperSolver is made visible.
function SuperSolver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SuperSolver (see VARARGIN)

% Choose default command line output for SuperSolver
handles.output = hObject;

%% Tabs Code
% Settings
TabFontSize = 8;
TabNames = {'Single Variable','Linear Systems','Nonlinear Systems'};
FigWidth = 0.35;

% Figure resize
set(handles.SimpleOptimizedTab,'Units','normalized')
pos = get(handles. SimpleOptimizedTab, 'Position');
set(handles. SimpleOptimizedTab, 'Position', [pos(1) pos(2) FigWidth pos(4)])

% Tabs Execution
handles = TabsFun(handles,TabFontSize,TabNames);

% Update handles structure
%handles.matrixInput = {}
handles.testvar = 0;
guidata(hObject, handles);

% UIWAIT makes SuperSolver wait for user response (see UIRESUME)
% uiwait(handles.SimpleOptimizedTab);

% --- TabsFun creates axes and text objects for tabs
function handles = TabsFun(handles,TabFontSize,TabNames)

% Set the colors indicating a selected/unselected tab
handles.selectedTabColor=get(handles.tab1Panel,'BackgroundColor');
handles.unselectedTabColor=handles.selectedTabColor-0.1;

% Create Tabs
TabsNumber = length(TabNames);
handles.TabsNumber = TabsNumber;
TabColor = handles.selectedTabColor;
for i = 1:TabsNumber
    n = num2str(i);
    
    % Get text objects position
    set(handles.(['tab',n,'text']),'Units','normalized')
    pos=get(handles.(['tab',n,'text']),'Position');
    
    % Create axes with callback function
    handles.(['a',n]) = axes('Units','normalized',...
        'Box','on',...
        'XTick',[],...
        'YTick',[],...
        'Color',TabColor,...
        'Position',[pos(1) pos(2) pos(3) pos(4)+0.01],...
        'Tag',n,...
        'ButtonDownFcn',[mfilename,'(''ClickOnTab'',gcbo,[],guidata(gcbo))']);
    
    % Create text with callback function
    handles.(['t',n]) = text('String',TabNames{i},...
        'Units','normalized',...
        'Position',[pos(3),pos(2)/2+pos(4)],...
        'HorizontalAlignment','left',...
        'VerticalAlignment','middle',...
        'Margin',0.001,...
        'FontSize',TabFontSize,...
        'Backgroundcolor',TabColor,...
        'Tag',n,...
        'ButtonDownFcn',[mfilename,'(''ClickOnTab'',gcbo,[],guidata(gcbo))']);
    
    TabColor = handles.unselectedTabColor;
end

% Manage panels (place them in the correct position and manage visibilities)
set(handles.tab1Panel,'Units','normalized')
pan1pos=get(handles.tab1Panel,'Position');
set(handles.tab1text,'Visible','off')
for i = 2:TabsNumber
    n = num2str(i);
    set(handles.(['tab',n,'Panel']),'Units','normalized')
    set(handles.(['tab',n,'Panel']),'Position',pan1pos)
    set(handles.(['tab',n,'Panel']),'Visible','off')
    set(handles.(['tab',n,'text']),'Visible','off')
end

% --- Callback function for clicking on tab
function ClickOnTab(hObject,~,handles)
m = str2double(get(hObject,'Tag'));

for i = 1:handles.TabsNumber;
    n = num2str(i);
    if i == m
        set(handles.(['a',n]),'Color',handles.selectedTabColor)
        set(handles.(['t',n]),'BackgroundColor',handles.selectedTabColor)
        set(handles.(['tab',n,'Panel']),'Visible','on')
    else
        set(handles.(['a',n]),'Color',handles.unselectedTabColor)
        set(handles.(['t',n]),'BackgroundColor',handles.unselectedTabColor)
        set(handles.(['tab',n,'Panel']),'Visible','off')
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = SuperSolver_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in singleVarListBox.
function singleVarListBox_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns singleVarListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from singleVarListBox


% --- Executes during object creation, after setting all properties.
function singleVarListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function singleVarToleranceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarToleranceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarToleranceEdit as text
%        str2double(get(hObject,'String')) returns contents of singleVarToleranceEdit as a double


% --- Executes during object creation, after setting all properties.
function singleVarToleranceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarToleranceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function singleVarIterationsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarIterationsEdit as text
%        str2double(get(hObject,'String')) returns contents of singleVarIterationsEdit as a double


% --- Executes during object creation, after setting all properties.
function singleVarIterationsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function singleVarx0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarx0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarx0Edit as text
%        str2double(get(hObject,'String')) returns contents of singleVarx0Edit as a double


% --- Executes during object creation, after setting all properties.
function singleVarx0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarx0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function singleVarFunctionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarFunctionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarFunctionEdit as text
%        str2double(get(hObject,'String')) returns contents of singleVarFunctionEdit as a double


% --- Executes during object creation, after setting all properties.
function singleVarFunctionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarFunctionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function newString = combineString(currentString, addedString)
newString = sprintf('%s\n%s',currentString, addedString);

% --- Executes on button press in singleVarSolveButton.
function singleVarSolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.singleVarOutputText, 'Max', 2);
method = get(handles.singleVarListBox, 'value');
syms infxn;
rawfxn = get(handles.singleVarFunctionEdit, 'string');
infxn = sym(rawfxn)
x0 = str2double(get(handles.singleVarx0Edit, 'string'));
Tol = str2double(get(handles.singleVarToleranceEdit, 'string'));
r = str2double(get(handles.singleVarRootEdit, 'string'));
iterations = str2double(get(handles.singleVarIterationsEdit, 'string'));
xList = zeros();
if isnan(r) || isempty(rawfxn) || isnan(x0) || isnan(Tol) || isnan(iterations)
    warning('off', 'symbolic:sym:sym:DeprecateExpressions');
    errorString = 'Missing Input';
    set(handles.singleVarOutputText, 'string', errorString);
    fprintf('Need Input')
else
    bandles = guidata(hObject)
    switch method
        case 1
            a = str2double(get(handles.bisectionAEdit, 'string'));
            b = str2double(get(handles.bisectionBEdit, 'string'));
            if isnan(a) || isnan(b)
                rangeErrorString='Missing range input';
                set(handles.singleVarOutputText, 'string', rangeErrorString);
            else
                aList = Bisection(infxn, a, b, r, Tol, bandles);
            end
        case 2
            aList = Single_Var_FixedPoint(infxn, x0, r, Tol, iterations, bandles);
        otherwise
            aList = Single_Var_Newtons(infxn, x0, r, Tol, iterations, bandles);
    end
    calcError(infxn, matlabFunction(diff(infxn)), r, aList(end), bandles);
end


function singleVarRootEdit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarRootEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarRootEdit as text
%        str2double(get(hObject,'String')) returns contents of singleVarRootEdit as a double


% --- Executes during object creation, after setting all properties.
function singleVarRootEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to singleVarRootEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==============================
% Single Variable Methods
%==============================

function xList = Single_Var_Newtons(infxn, x0, r, Tol, iterations, handles)

xList = zeros();                        % List of x values calculated for graphing and tables
xList(1) = x0;                          % Set the first element of the list to the initial guess
iterativeErrorList = zeros();           % List of ei
errorFlag = false;                      % divergence flag
f = matlabFunction(infxn);              % Convert the symbolic function to a function handle
fprime = matlabFunction(diff(infxn));   % First derivative of f
count = zeros();

% For analysis
pastError = 0;
dx = 1;

if(fprime(r) == 0)
    fprintf('f''(r) = 0 so Newton''s Method will be linearly convergent\n');
    linearString = '\n f''(r) = 0 so Newton''s Method will be linearly convergent.';
    set(handles.singleVarOutputText, 'string', linearString);
else
    fprintf('The method will converge quadratically.');
    quadConString = 'The method will converge quadratically.';
    set(handles.singleVarOutputText, 'string', quadConString);
end

for i = 1:iterations
    disp('iteration: ');
    disp(i);
    count(i) = i;
    tic;
    fofx = f(xList(i));
    fpofx = fprime(xList(i));
    
    if(fpofx == 0 || abs(fpofx) == Inf)
        errorString = ['The derivative of the function at ', num2str(xList(i)), ' is 0, try another initial guess'];
        set(handles.singleVarOutputText, 'string', errorString);
        errorFlag = true;
        break;
    end
    
    xList(i+1) = xList(i) - fofx/fpofx ;             % Gets the next value of x
    ei = abs(xList(i) - r);                    % Gets forward error of current iteration
    iterativeErrorList(i) = ei/(pastError)^2;  % Error in relation to previous iteration
    
    % Checks if the difference in x values has converged
    if (dx <= Tol || abs(fofx) <= Tol)
        break;
    end
    % Checks if error still exists
    if iterativeErrorList(i) == Inf || iterativeErrorList(i) > 1  % Handles case where result is invalid
        iterativeErrorList(i) = NaN;
    end
    
    % Prep for next iteration of loop
    pastError = ei;
end
count(i+1) = i+1;
figure
plot(count, xList);
xlabel('Iteration');
ylabel('f(x)');
title('Newton Method of Single Variable');
time = toc;
timeString = ['Elapsed time = ', num2str(time)];
set(handles.singleVarOutputText, 'string', timeString);
fprintf('root = %12.8f\n', xList(i));
% Fixed Point Iteration
% Receives:
% infxn = function in x
% x0 = initial guess
% iterations = number of iterations
% Returns list of calculated x values

function xList = Single_Var_FixedPoint(infxn, x0, r, Tol, iterations, handles)
xList = zeros;          % Have x be a list for creating graphs
count = zeros;          % Iteration array for creating graphs

f = matlabFunction(infxn);
fprime = matlabFunction(diff(infxn));
% TODO check for convergence
fprimeofr = fprime(r);
if abs(fprimeofr) > 1
    fprintf('May not converge')
    set(handles.singleVarOutputText, 'string', 'FPI may not converge')
elseif fprimeofr == 0
    fprintf('Quadratically convergent')
    set(handles.singleVarOutputText, 'string', 'Since g''(r) = 0, FPI will be quadratically convergent')
else
    fprintf('Linearly convergent with rate %s', num2str(abs(fprimeofr)))
    set(handles.singleVarOutputText, 'string', 'Since |g''(r)| < 1, FPI will be quadratically convergent')
end
xList(1) = x0;
% Runs until root approximated or number of iterations reached
for i = 1:iterations
    count(i) = i;
    tic;
    xList(i+1) = f(xList(i));           % Get next x from result of current x
    dx = abs(xList(i+1) - xList(i));    % Tracks amount result changed
    if abs(dx) < Tol
        break;
    end
end
time = toc;
count(i+1) = i+1;

xa = xList(i);

figure
plot(count, xList);
xlabel('Iteration');
ylabel('g(x)');
title('Fixed Point Iteration');
timeString = ['Elapsed time = ', num2str(time)];
set(handles.singleVarOutputText, 'string', timeString);
fprintf('Approximate root = %8f\n', xa);

%Program 1.1 Bisection Method
%Computes approximate solution of f(x)=0
%Input: function handle f; a,b such that f(a)*f(b)<0,
% and tolerance tol
%Output: Approximate solution xc
function cList=Bisection(infxn,a,b,r,tol, handles)
syms x;
cList = zeros;
count = zeros;          % Iteration array for creating graphs


f = matlabFunction(infxn);
tic;
if sign(f(a))*sign(f(b)) >= 0
    errorString = 'f(a)f(b)<0 not satisfied!'; %ceases execution
    set(handles.singleVarOutputText, 'string', errorString);
else
    fa=f(a);
    i = 1;
    while (b-a)/2>tol
        cList(i)=(a+b)/2;
        fc=f(cList(i));
        if fc == 0 %c is a solution, done
            break
        end
        if sign(fc)*sign(fa)<0 %a and c make the new interval
            b=cList(i);%fb=fc;
        else %c and b make the new interval
            a=cList(i);fa=fc;
        end
        i = i + 1;
        count(i) = i;
    end
    cList(i)=(a+b)/2; %new midpoint is best estimate

end
time = toc;

figure
plot(count, cList);
xlabel('Iteration');
ylabel('c');
title('Bisection');
timeString = ['Elapsed time = ', num2str(time)];
set(handles.singleVarOutputText, 'string', timeString);

% function to calculate forward error, backward error, and error
% magnification of methods used to solve single variable equations
% Pass in original syms function, its derivative, the real root, and the
% approximate root
function calcError(func, deriv, r, xa, handles)
% might have to remove error magnification or find alternative for gpow
%gPow = feval(symengine, 'degree', func);    % Gets g(x) from highest degree of the equations, func has to be symbolic for degree() to work
%gofr = r^gPow;                              % g(r) will just be r^degree
%magError = abs(gofr/(r*feval(deriv, r)));   % Calculating the magnitude of error from equation on page 49
forwardErr = abs(r - xa);
backwardErr = double(abs(subs(func,xa)));
rootM = getRootMultiplicity(func, r);
approxStr = 'Approximate Root = ';
realString = 'Real Root = ';
mString = ' has multiplicity = ';
fString = 'Forward Error = ';
bString = 'Backward Error = ';
currString = get(handles.singleVarOutputText, 'string');
newString = sprintf('%s%s\n%s%s%s%s\n%s%s\n%s%s\n', approxStr, num2str(xa), realString, num2str(r), mString, num2str(rootM), fString, num2str(forwardErr), bString, num2str(backwardErr));
statusString = combineString(currString, newString);
set(handles.singleVarOutputText, 'string', statusString);
% Remove and return the error magnifcation when integrating code
fprintf('Real Root = %1.8f \nApproximate Root = %12.8f\nForward Error = %12.8f \nBackward Error = %12.8f\nError Magnification = %12.8f\n', r, xa, forwardErr, backwardErr,1);

function multiplicity = getRootMultiplicity(fnx, x)
    multiplicity = 0;
    while(subs(fnx, x) == 0)        % Runs when fn(0) = 0
        multiplicity = multiplicity + 1;
        fnx = diff(fnx);            % Get the next derivative
    end

function bisectionAEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bisectionAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bisectionAEdit as text
%        str2double(get(hObject,'String')) returns contents of bisectionAEdit as a double


% --- Executes during object creation, after setting all properties.
function bisectionAEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bisectionAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bisectionBEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bisectionBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bisectionBEdit as text
%        str2double(get(hObject,'String')) returns contents of bisectionBEdit as a double


% --- Executes during object creation, after setting all properties.
function bisectionBEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bisectionBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in linearSysListBox.
function linearSysListBox_Callback(hObject, eventdata, handles)

% hObject    handle to linearSysListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns linearSysListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from linearSysListBox


% --- Executes during object creation, after setting all properties.
function linearSysListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linearSysListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function linearSysNumOfVarEdit_Callback(hObject, eventdata, handles)
% hObject    handle to linearSysNumOfVarEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linearSysNumOfVarEdit as text
%        str2double(get(hObject,'String')) returns contents of linearSysNumOfVarEdit as a double


% --- Executes during object creation, after setting all properties.
function linearSysNumOfVarEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linearSysNumOfVarEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function linearSysIterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to linearSysIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linearSysIterEdit as text
%        str2double(get(hObject,'String')) returns contents of linearSysIterEdit as a double


% --- Executes during object creation, after setting all properties.
function linearSysIterEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linearSysIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function linearSysTolEdit_Callback(hObject, eventdata, handles)
% hObject    handle to linearSysTolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linearSysTolEdit as text
%        str2double(get(hObject,'String')) returns contents of linearSysTolEdit as a double


% --- Executes during object creation, after setting all properties.
function linearSysTolEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linearSysTolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function omegaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to omegaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of omegaEdit as text
%        str2double(get(hObject,'String')) returns contents of omegaEdit as a double


% --- Executes during object creation, after setting all properties.
function omegaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omegaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in linearSysMatrixButton.
function linearSysMatrixButton_Callback(~, ~, handles)
global tableInput
rows = str2double(get(handles.linearSysNumOfVarEdit, 'string'));

if isnan(rows) || rows < 1 || rows == Inf
    errorString = 'Invalid variable amount';
    set(handles.lSysOutputText, 'string', errorString);
else
    dataInput = cell(rows, rows + 1);
    figure,
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1:2) [660 120]])
    %Input table's creation
    tableInput = uitable('ColumnWidth',{50},...
        'Position',[30 20 600 80], ...
        'data',dataInput, ...
        'ColumnEditable',true);
end
%gen data because so i do not have to enter it in manually for the
%example
%handles.matrixInput = get(tableInput,'data')
% hObject    handle to linearSysMatrixButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in linearSysSolveButton.
function linearSysSolveButton_Callback(hObject, ~, handles)
linHandles = guidata(hObject);
global tableInput
global initialMatrixInput

retrievedData = get(tableInput, 'data');
augA = str2double(retrievedData);
method = get(handles.linearSysListBox, 'value');
Tol = str2double(get(handles.linearSysTolEdit, 'string'));
omega = str2double(get(handles.omegaEdit, 'string'));
iterations = str2double(get(handles.linearSysIterEdit, 'string'));
set(handles.lSysOutputText, 'Max', 2);

if isnan(Tol) || isnan(iterations)
    fprintf('Need Input')
else
    switch method
        case 1
            sol = Gauss_Elim(augA, linHandles);
        case 2
            sol = LU_Decomposition(augA, linHandles)
        case 3
            retrievedGuess = get(initialMatrixInput, 'data')
            P = str2double(retrievedGuess)
            sol = Jacobi(augA, P, linHandles)
        otherwise
            if isnan(omega)
                errStr = 'need input';
                set(handles.lSysOutputText, 'string', errStr);
            else
                retrievedGuess = get(initialMatrixInput, 'data')
                x0 = str2double(retrievedGuess);
                sol = SOR(augA, x0, omega, iterations, Tol, linHandles)
            end
    end
end


% hObject    handle to linearSysSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in linearSysGuessButton.
function linearSysGuessButton_Callback(hObject, eventdata, handles)
set(handles.lSysOutputText, 'Max', 2);
vars = str2double(get(handles.linearSysNumOfVarEdit, 'string'))
dataInput = cell(1, vars)

global initialMatrixInput
figure,
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) [660 120]])
%Input table's creation
initialMatrixInput = uitable('ColumnWidth',{70},...
    'Position',[30 20 600 80], ...
    'data',dataInput, ...
    'ColumnEditable',true);
% hObject    handle to linearSysGuessButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%==============================
% Linear Systems Methods
%==============================

% Naive Gaussian Elimination
% x,y,z variables specification
% # of rows
% # of cols
% row input csv values?
function solutions = Gauss_Elim(augA, handles)
opCount = 0;
%row = 3;
%col = 3;

% CONDITION NUMBER
% Conditioning is a property of the matrix
% Error Magnification factors of the magnitude cond(matrix) are possible
% Matlab default precision is double


% TODO incorporate user input
%A = [1, 2, -1; 2, 1, -2; -3, 1, 1]
%ans = [3; 3; -6];
solutions = zeros;
n = size(augA, 1);
row = size(augA, 1);
col = row + 1;

A = augA;
A(:,col) = [];
norminf = norm(A, inf);
norminf_inv = norm(inv(A), inf);

cond_num = norminf * norminf_inv;

iso_exp = floor(log10(cond_num*10));
fprintf('Error Magnification factors of the magnitude %d are possible.\n', iso_exp);
fprintf('Since Matlab defaults to double precision this means that \n');
fprintf('16 - %d = %d correct digits in the solution.\n', iso_exp, 16-iso_exp);

condStr1 = ['Error Magnification factors of the magnitude ', num2str(iso_exp), ' are possible'];
condStr2 = 'Since Matlab defaults to double precision this means that';
condStr3 = ['16 - ', num2str(iso_exp),' = ', num2str(16-iso_exp), ' correct digits in the solution.'];

condStr = sprintf('%s\n%s\n%s\n', condStr1, condStr2, condStr3);

set(handles.lSysOutputText, 'string', condStr);

tic;

for j = 1:n-1 % n-1 = num of rows - the 1st row
    % check for division by zero
    if augA(j,j) == 0
        opCount = opCount + 1;
        error('zero pivot encountered');
    end
    % eliminate col j to put zero in each location below diag
    % ex. a(j+1, j), a(j+2, j),...a(n,j)
    for i = j+1:n
        % row multiplier
        multi = augA(i, j) / augA(j, j);
        opCount = opCount + 1;

        % subtract multiplier * the row from 
        for index = 1:col
            augA(i, index) = (augA(i, index) - (multi * augA(j, index)))
            opCount = opCount + 1;
        end
    end
end

% Backsolving
for q = n:-1 : 1
    solutions(q) = augA(q, row + 1);
    for u = q+1:n
        opCount = opCount + 1;
        solutions(q) = solutions(q) - augA(q,u)*solutions(u);
    end
    solutions(q) = solutions(q)/augA(q,q);
    fprintf('\n');
end
    
for x = 1:n
    fprintf('x%d = %f\n', x, solutions(x));
end
fprintf('\n');
time = toc;
fprintf('\n');

opsString = ['Number of Operations = ', num2str(opCount), ' seconds'];
timeString = ['Elapsed Time = ', num2str(time)];
newString = combineString(condStr, opsString);
newString = combineString(newString, timeString);
set(handles.lSysOutputText, 'string', newString);
%fprintf('Number of Operations = %d\n\n', opCount);

% This function performs the LU factorization of a symmetric matrix.
% The function checks for singularity by checking for a zero pivot.
% The function also does a partial pivoting protocol to prevent swamping.
% This function also creates a permutation matrix to help with backsolving.
%   Input: Symmetric matrix and solution vector 
%   Output: Lower Triangular factor
%           Upper Triangular factor
%           Permutation matrix
%           Solution Vector x

% input from user for now. may need to update.
function solution = LU_Decomposition(augA, handles)
%augA = [1 0 1; 0 1 1, 
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
%A = input('enter an initial matrix: \n');
b = augA(:,columns);
opCount = 0;


% % % matrix A
% A=[-1, 0, 1;
%    2, 1, 1;
%   -1, 2, 0]
 %vector b
% b = [-2; 17; 3]

% size of n x n
[n, n] = size(A);

% init U to A
U = A;

% init L to identity
L = eye(n);

% identity matrix for pivoting
P = eye(n);

tic;
for k = 1:n
    % create matrix pivot
    %  maximum absolute value from the k thru n rows, k column
    [pivot i] = max(abs(U(k:n,k)));
    i = i + k-1;
    opCount = opCount + 1;
    % if i != k
    if i ~= k 
        % interchange rows i and k in U
        % U(k, :) means (row k, all columns in row k)
        opCount = opCount + 3;
        temp = U(k, :);
        U(k, :) = U(i, :);
        U(i, :) = temp;
        % interchange rows i and k in P
        opCount = opCount + 3;
        temp = P(k, :);
        P(k, :) = P(i, :);
        P(i, :) = temp;
        
        % Lower triangular array
        if k >= 2
            % (k , 1:k-1) means entire row k and columns 1 thru k-1
            temp = L(k, 1:k-1);
            L(k, 1:k-1) = L(i, 1:k-1);
            L(i, 1:k-1) = temp;
            opCount = opCount + 3;
        end %endif
    end %end outer if ~=
    % get factors L and U
    for j = k+1:n
        % set multiplier in L
        L(j, k) = U(j, k) / U(k, k); 
        opCount = opCount + 1;
        % for all in row j = for all in row j - multiplier*for all in row k
        U(j, :) = U(j, :) - L(j, k) * U(k, :);
        opCount = opCount + 1;
    end % end gauss elim
end

%solve for c
c = [];
inv(L)
Pb = P*b;
c = L\Pb;
opCount = opCount + n^2;

%solve for x
solution = [];
solution = U\c;
opCount = opCount + n^2;

time = toc;

opsString = ['Number of Operations = ', num2str(opCount), ' seconds'];
timeString = ['Elapsed Time = ', num2str(time)];
newString = combineString(opsString, timeString);
set(handles.lSysOutputText, 'string', newString);

fprintf('Number of Operations = %d\n', opCount);
fprintf('Lower Triangular = \n');
L
fprintf('Upper Triangular = \n');
U
fprintf('Permutation Matrix = \n');
P
fprintf('Solution vector = \n');
solution

function X = Jacobi(augA, P, handles)
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
%A = input('enter an initial matrix: \n');
b = augA(:,columns);
opCount = 0;
convergenceStr = '';

X = zeros();
iterations = 10;
N = length(b);
Tol = 0.00000001;
diagonalDominant = true;

% Check if A is strictly diagonally dominant
% For each row in matrix A
for r=1:N
    rowSum = sumabs(A(r,:)) - abs(A(r, r)); % Sum of the entire row minus the diagonal
    % Check if diagonal is strictly greater than the row sum
    if abs(A(r,r)) < rowSum
        diagonalDominant = false;           % If not, note ite
        break;
    end
end

if diagonalDominant == false
    fprintf('Matrix A is not strictly diagonal dominant and may not converge.\n');
    convergenceStr = 'Matrix A is not strictly diagonal dominant and may not converge.';
    set(handles.lSysOutputText, 'string', convergenceStr);
end

tic;

for k=1:iterations
    for j=1:N
        X(j)=(b(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);     % slicing using coefficients to solve for uk and vk
    end
    err=abs(norm(X'-P));
    relerr=err/(norm(X)+eps);
    P=X';
    if (err<Tol)||(relerr<Tol)
        break
    end
end

time = toc;

opsString = ['Number of Operations = ', num2str(opCount)];
timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
newString = combineString(opsString, timeString);
set(handles.lSysOutputText, 'string', newString);

function[x, error, iter, flag]  = SOR(augA, x, w, max_it, tol, handles)
rows = size(augA,1);
columns = rows + 1;
A = augA;
A(:,columns) = [];
%A = input('enter an initial matrix: \n');
b = augA(:,columns);

opCount = 0;

convStr = '';
% input initialization
% A = [3 1 -1; 
%      2 4 1; 
%     -1 2 5];
% 
% [n, n] = size(A);
% 
% xold = [0; 0; 0];
% b = [4; 1; 1];
% w = 1.25;
% max_it = 2;
% tol = 0.00000001;


x = x';
diagonalDominant = true;

n = length(A);

if (w == 1)
    fprintf('The relaxation scalar omega = 1. The method used is now Gauss-Seidel')
    GSString = 'The relaxation scalar omega = 1. The method used is now Gauss-Seidel';
    set(handles.lSysOutputText, 'string', GSString);
end

tic;
% check for diagonally dominant convergence guarantee
for r=1:n 
    % Sum of the entire row minus the diagonal
    rowSum = sumabs(A(r,:)) - abs(A(r, r)); 
    opCount = opCount + 1;
    % Check if diagonal is strictly greater than the row sum
    if (abs(A(r,r)) < rowSum)
    % If not, note it
        diagonalDominant = false;           
        break;
    end
end

if (diagonalDominant == false)
% Let user know that convergence is not guaranteed
    fprintf('Matrix A is not strictly diagonal dominant and may not converge.\n');
    convStr = 'Matrix A is not strictly diagonal dominant and may not converge';
    newStr = combineString(GSString, convStr);
    set(handles.lSysOutputText, 'string', newStr);
end

if size(b) ~= size(x)
    fprintf('The given approximation vector does not match the x vector size\n');
    errStr = 'The given approximation vector does not match the x vector size';
    set(handles.lSysOutputText, 'string', errStr);
end


flag = 0;    
count = 1;

% matrix splitting 
% TODO: opCount will be affected by these operations. Need to add.
D = diag(diag(A));
L = tril(A-D);
U = triu(A-D);
%M = D + w*L;
%N = (1 - w)*D - w*U;
%G = M\N;
G = (inv(D+w*L))*(((1-w)*D-w*U)*x +w*b);
opCount = opCount + 1;
%fprintf('Iteration 1: ', G);
datasave = [];
% begin iteration
for iter = 1:max_it                         
        xnew = G;
        RelForError = (norm(xnew-x))/(norm(xnew));
        opCount = opCount + 1;
        % update approximation
        while (RelForError > tol)
            x = xnew;
            opCount = opCount + 1;
            G = (inv(D+w*L))*(((1-w)*D-w*U)*x +w*b)
            xnew = G;
            opCount = opCount + 1;
            RelForError = (norm(xnew-x))/(norm(xnew));
            if (RelForError <= tol)
                break
            end
            count = count+1;
            x = [x, xnew];
            datasave = [datasave; count, RelForError, flag];
        end
end

b = b / w; % vector b
if (RelForError > tol) 
   flag = 1;
   fprintf('Did not converge')
   convStr = 'Did not converge';
   set(handles.lSysOutputText, 'string', convStr);
end

% for function return
x = xnew;
error = RelForError;
iter = count;

fprintf('Number of operations: %d\n', opCount);
time = toc;

opsStr = ['Number of Operations = ', num2str(opCount)];
timeString = ['Elapsed Time = ', num2str(time), ' seconds'];
newStr = combine(convStr, opsStr);
newStr = combineString(newStr, timeString);
set(handles.lSysOutputText, 'string', newStr);

fprintf('\n');

fprintf('  iteration    error    flag\n')

disp(datasave)
fprintf(' x final\n')
disp(xnew)

% -------------------------
% --- NONLINEAR SYSTEMS ---
% -------------------------

% --- Executes on selection change in nonLinearSysListBox.
function nonLinearSysListBox_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nonLinearSysListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nonLinearSysListBox


% --- Executes during object creation, after setting all properties.
function nonLinearSysListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSysInitialGuessEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysInitialGuessEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysInitialGuessEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysInitialGuessEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysInitialGuessEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysInitialGuessEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSysEqnListEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysEqnListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysEqnListEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysEqnListEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysEqnListEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysEqnListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSySIterationsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSySIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSySIterationsEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSySIterationsEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSySIterationsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSySIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BroydenInitialMatrix.
function BroydenInitialMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to BroydenInitialMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BroydenInitialMatrix
% Hint: popupmenu controls usually have a white background on Windows.



function nonLinearSysVarsListEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysVarsListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysVarsListEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysVarsListEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysVarsListEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysVarsListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSysBroydenEqnListEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenEqnListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysBroydenEqnListEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysBroydenEqnListEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysBroydenEqnListEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenEqnListEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSysBroydenGuessEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenGuessEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysBroydenGuessEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysBroydenGuessEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysBroydenGuessEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenGuessEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonLinearSysBroydenIterationsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysBroydenIterationsEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysBroydenIterationsEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysBroydenIterationsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenIterationsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nonLinearMultiVarSolveButton.
function nonLinearMultiVarSolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearMultiVarSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonLinearMultiVarSolveButton

%get variable list
variables = get(handles.nonLinearSysVarsListEdit, 'string');
ca = {};

%remove white spaces and convert string to char array
variables(variables==' ') = [];
variables = char(variables); 

%load each variable into a cell
for i=1:length(variables)
    ca{i} = variables(i); 
end
celldisp(ca);

%convert each char to sym type
vars = cell2sym(ca);
for i=1:length(variables)
     syms(variables(i))
end

%GET initial guess
x = str2num(get(handles.nonLinearSysInitialGuessEdit, 'string'));
disp(x);

%GET system of equations
eqn_str = get(handles.nonLinearSysEqnListEdit, 'string');
eqns = sym(eqn_str);
disp(eqns);

%GET number of iterations
number_of_iterations = str2double(get(handles.nonLinearSySIterationsEdit, 'string'));
disp(number_of_iterations);

%pass parameters to MultiVarNewton
Multi_Var_Newton_Method(vars, x, eqns, number_of_iterations);

function  Multi_Var_Newton_Method(vars, x, eqns, number_of_iterations)
%create a Jacobian matrix
DF = jacobian(eqns, vars);
disp(DF);
x_values = zeros(1,100);
y_values = zeros(1,100);
t = zeros(1,100);
%begin iteration steps for calculating the solution of the system
for i=1:number_of_iterations
    disp('iteration: ');
    disp(i);
    %solve for the solution set, s, to plug into later
    tic;
    a = zeros(length(eqns),1);
  
    for j=1:length(eqns)
        answer = subs(eqns(j), vars, x);
        a(j) = single(answer);
    end %end of solution set loop

   
    %solve for the constants in the Jacobian matrix
    %find values of Jacobian with starting point
    sol_matrix = subs(DF, vars, x);

    %a = array containing the solution variables [s1, s2, s3...],
    %sol_matrix = solution matrix to solve for the solution set (s1, s2, s3...)
    %solve the linear system of equations this will output the solution set
    %{'a', 'b'}
    %[1,1]

    %[(a+b), (a+(-b)), (a+b^2)]
    %
    
    sol_set = linsolve(sol_matrix, a);
    if(isinf(sol_set))
        warning('The system of equations does not converge')
        return
    end
    sol_set = double(sol_set);
    
    %reshape to a 1D array for easy subtraction
    sol_set = reshape(sol_set, [1, numel(a)]);


    %solve xk = x(k-1) + s --> xk
    x_values(i) = x(1);
    y_values(i) = x(2);
    prev_x = x;
    disp(vpa(x,10))
    x = x - sol_set;
    if((round(sum(prev_x - x), 16)) == 0)
        disp('solution found at x = ');
        disp(vpa(x))
        figure
        subplot(2,1,1)       % add first plot in 2 x 1 grid
        plot(x_values,y_values)
        title('Convergence or Divergence')

        subplot(2,1,2)       % add second plot in 2 x 1 grid
        plot(t)       % plot using + markers
        title('TIme Complexity')
        drawnow()
        return
    end
    diverge = 0; 
    if(i <= 10)
        distance(i) = round(sum(prev_x - x));
        if(i == 10)
            if(distance(i) > distance(i-9))
                diverge = 1;
            end
        end
    end
    
    
    %show the solution at the end of each step[[2*u^2 + v^2 + 3*w^2 + 6*w - 4*u + 2],[3*u^2 - 12*u + v^2 + 3*w^2 + 8], [u^2 + v^2 - 2*v + 2*w^2 - 5]]
   t(i) = toc;
end %repeat k times end of iteration loop

if(diverge == 1)
    disp('the solution appears to be diverging')
    
    
end
% plot(x_values, y_values)
% drawnow()
% plot(t)
% drawnow()
figure
subplot(2,1,1)       % add first plot in 2 x 1 grid
plot(x_values,y_values)
title('Convergence or Divergence')

subplot(2,1,2)       % add second plot in 2 x 1 grid
plot(t)       % plot using + markers
title('Time Complexity')
    


% --- Executes on button press in nonLinearSysBroyden1SolveButton.
function nonLinearSysBroyden1SolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroyden1SolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonLinearSysBroyden1SolveButton
x_i = str2num(get(handles.nonLinearSysBroydenGuessEdit, 'string'));
disp(x_i);
eqns = get(handles.nonLinearSysBroydenEqnListEdit, 'string');
eqns = strsplit(eqns, ';');
for i=1:length(eqns)
    eqns{i} = strtrim(eqns{i});
    eqns{i} = str2func(eqns{i});
end
disp(eqns);
A = str2num(get(handles.nonLinearSysBroydenMatrixEdit, 'string'));
disp(A);
number_of_iterations = str2double(get(handles.nonLinearSysBroydenIterationsEdit, 'string'));
disp(number_of_iterations);
Broyden_1_Method(x_i, eqns, A, number_of_iterations);

function Broyden_1_Method(x_i, eqns, A, number_of_iterations)
%evaluate the functions at x
%{(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
x = num2cell(x_i);
y1 = cellfun(@(t) t(x{:}), eqns);
x1 = x_i;
x1 = transpose(x1);
y1 = transpose(y1);
x_values = zeros(1,100);
y_values = zeros(1,100);
t = zeros(1,100);
%begin iteration steps
for i=1:number_of_iterations
    x = x1;
    x_values(i) = x(1);
    y_values(i) = x(2);
    y = y1;
    tic;
    x1 = x - A\y;
    disp(x1);
    x_val = num2cell(x1);
    
    y1 = cellfun(@(t) t(x_val{:}), eqns);
    if(round(y1, 10) == 0)
        disp('solution found at x = ')
        disp(vpa(x1,10))
        disp('correct to 10 decimals digits')
        plot(t)
        title('Time Complexity Graph')
        xlabel('Iterations')
        ylabel('Time to Compute (s)')
        return
    end
    y1 = transpose(y1);
    
    %calculate new matrix A
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = A*deltaX;
    p2 = deltaY - p1;
    p3 = p2 * transpose(deltaX);
    p4 = transpose(deltaX)*deltaX;
    A = A + p3/p4;
    t(i) = toc;
end
figure, plot(t)



function nonLinearSysBroydenMatrixEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenMatrixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonLinearSysBroydenMatrixEdit as text
%        str2double(get(hObject,'String')) returns contents of nonLinearSysBroydenMatrixEdit as a double


% --- Executes during object creation, after setting all properties.
function nonLinearSysBroydenMatrixEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroydenMatrixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nonLinearSysBroyden2SolveButton.
function nonLinearSysBroyden2SolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysBroyden2SolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonLinearSysBroyden2SolveButton
x_i = str2num(get(handles.nonLinearSysBroydenGuessEdit, 'string'));
disp(x_i);
eqns = get(handles.nonLinearSysBroydenEqnListEdit, 'string');
eqns = strsplit(eqns, ';');
for i=1:length(eqns)
    eqns{i} = strtrim(eqns{i});
    eqns{i} = str2func(eqns{i});
end
disp(eqns);
A = str2num(get(handles.nonLinearSysBroydenMatrixEdit, 'string'));
disp(A);
number_of_iterations = str2double(get(handles.nonLinearSysBroydenIterationsEdit, 'string'));
disp(number_of_iterations);
Broyden_2_Method(x_i, eqns, A, number_of_iterations);

function Broyden_2_Method(x_i, eqns, B, number_of_iterations)
%evaluate the functions at x
%{(@(u,v)u^2+v^2-1) ,   (@(u,v)(u-1)^2+v^2-1)}
x = num2cell(x_i);
y1 = cellfun(@(t) t(x{:}), eqns);
x1 = x_i;
x1 = transpose(x1);
y1 = transpose(y1);
t = zeros(1,100);
%begin iteration steps
for i=1:number_of_iterations
    x = x1;
    y = y1;
    tic;
    x1 = x - B*y;
    x_val = num2cell(x1);
    display(x1);
    y1 = cellfun(@(t) t(x_val{:}), eqns);
    if(round(y1, 10) == 0)
        disp('solution found at x = ')
        disp(vpa(x1,10))
        disp('correct to 10 decimals digits')
        plot(t)
        title('Time Complexity Graph')
        xlabel('Iterations')
        ylabel('Time to Compute (s)')
        drawnow()
        return
    end
    y1 = transpose(y1);
    
    %calculate new matrix B
    deltaY = y1 - y;
    deltaX = x1 - x;
    p1 = B*deltaY;
    p2 = deltaX - p1;
    p3 = p2 * (transpose(deltaX)*B);
    p4 = transpose(deltaX)*B*deltaY;
    B = B + p3/p4;
    t(i) = toc;
end
figure, plot(t)
%>>>>>>> origin/master



function lSysOutputText_Callback(hObject, eventdata, handles)
% hObject    handle to lSysOutputText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lSysOutputText as text
%        str2double(get(hObject,'String')) returns contents of lSysOutputText as a double


% --- Executes during object creation, after setting all properties.
function lSysOutputText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lSysOutputText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
