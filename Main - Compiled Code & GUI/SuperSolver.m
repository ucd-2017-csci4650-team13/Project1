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


% --- Executes on button press in singleVarSolveButton.
function singleVarSolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('Button runs')
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
    fprintf('Need Input')
else
    switch method
        case 1
            a = str2double(get(handles.bisectionAEdit, 'string'))
            b = str2double(get(handles.bisectionBEdit, 'string'))
            if isnan(a) || isnan(b)
                error('Need range input')
            end
            xList = Bisection(infxn, a, b, r, Tol);
        case 2
            xList = Single_Var_FixedPoint(infxn, x0, r, Tol, iterations);
        otherwise
            xList = Single_Var_Newtons(infxn, x0, r, Tol, iterations);
    end
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

function xList = Single_Var_Newtons(infxn, x0, r, Tol, iterations)
%Tol = 0.00000001;
%syms x;
xList = zeros();                        % List of x values calculated for graphing and tables
xList(1) = x0;                          % Set the first element of the list to the initial guess
iterativeErrorList = zeros();           % List of ei
%infxn = sym(infxn)
f = matlabFunction(infxn);              % Convert the symbolic function to a function handle
fprime = matlabFunction(diff(infxn));   % First derivative of f

% For analysis
pastError = 0;
dx = 1;

if(fprime(r) == 0)
    fprintf('f''(r) = 0 so Newton''s Method will be linearly convergent\n');
else
    fprintf('The method will converge quadratically.\n');
end
for i = 1:iterations
    fofx = f(xList(i));
    fpofx = fprime(xList(i));
    
    if(fpofx == 0 || abs(fpofx) == Inf)
        error('The derivative of the function at %12.8f is 0, try another initial guess', xList);
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

fprintf('root = %12.8f\n', xList(i));
calcError(infxn, fprime, r, xList(i));
fprintf('root %12.8f has m = %8f\n', r, getRootMultiplicity(infxn, r));

% Fixed Point Iteration
% Receives:
% infxn = function in x
% x0 = initial guess
% iterations = number of iterations
% Returns list of calculated x values

function xList = Single_Var_FixedPoint(infxn, x0, r, Tol, iterations)
%Tol = 0.00000001;       % Stopping criteria
xList = zeros;          % Have x be a list for creating graphs

f = matlabFunction(infxn);
fprime = matlabFunction(diff(infxn));
% TODO check for convergence

xList(1) = x0;
% Runs until root approximated or number of iterations reached
for i = 1:iterations
    xList(i+1) = f(xList(i));           % Get next x from result of current x
    dx = abs(xList(i+1) - xList(i));    % Tracks amount result changed
    if abs(dx) < Tol
        break;
    end
end

xa = xList(i);
calcError(infxn, fprime, r, xList(i));

fprintf('Approximate root = %8f\n', xa);

%Program 1.1 Bisection Method
%Computes approximate solution of f(x)=0
%Input: function handle f; a,b such that f(a)*f(b)<0,
% and tolerance tol
%Output: Approximate solution xc
function xc=Bisection(infxn,a,b,r,tol)
syms x;
f = matlabFunction(infxn);
if sign(f(a))*sign(f(b)) >= 0
    error('f(a)f(b)<0 not satisfied!') %ceases execution
end
fa=f(a);
fb=f(b);
i = 1;
cList = zeros();
while (b-a)/2>tol
    cList(i)=(a+b)/2;
    fc=f(cList(i));
    if fc == 0 %c is a solution, done
        break
    end
    if sign(fc)*sign(fa)<0 %a and c make the new interval
        b=cList(i);fb=fc;
    else %c and b make the new interval
        a=cList(i);fa=fc;
    end
    i = i + 1;
end
xc=(a+b)/2; %new midpoint is best estimate
calcError(infxn, matlabFunction(diff(infxn)), r, xc)

% function to calculate forward error, backward error, and error
% magnification of methods used to solve single variable equations
% Pass in original syms function, its derivative, the real root, and the
% approximate root
function calcError(func, deriv, r, xa)
% might have to remove error magnification or find alternative for gpow
%gPow = feval(symengine, 'degree', func);    % Gets g(x) from highest degree of the equations, func has to be symbolic for degree() to work
%gofr = r^gPow;                              % g(r) will just be r^degree
%magError = abs(gofr/(r*feval(deriv, r)));   % Calculating the magnitude of error from equation on page 49
forwardErr = abs(r - xa);
backwardErr = abs(subs(func,xa));
% Remove and return the error magnifcation when integrating code
%fprintf('Real Root = %1.8f \nApproximate Root = %12.8f\nForward Error = %12.8f \nBackward Error = %12.8f\nError Magnification = %12.8f\n', r, xa, forwardErr, backwardErr, magError);



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
function linearSysMatrixButton_Callback(hObject, eventdata, handles)
rows = str2double(get(handles.linearSysNumOfVarEdit, 'string'))
dataInput = cell(rows, rows + 1)
%dataInput = cell(2,3)
global retrievedData
global tableInput
figure,
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2) [660 120]])
%Input table's creation
tableInput = uitable('ColumnWidth',{70},...
    'Position',[30 20 600 80], ...
    'data',dataInput, ...
    'ColumnEditable',true);
%gen data because so i do not have to enter it in manually for the
%example
retrievedData = get(tableInput,'data')
%handles.matrixInput = get(tableInput,'data')
% hObject    handle to linearSysMatrixButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in linearSysSolveButton.
function linearSysSolveButton_Callback(hObject, eventdata, handles)
global retrievedData
global tableInput
retrievedData = get(tableInput, 'data')
g= cell2mat(retrievedData)
k = handles.testvar
i = g(1,1)
%g = get(handles.matrixInput, 'data')
% hObject    handle to linearSysSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
