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
FigWidth = 0.5;

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
global rows 
rows = str2double(get(handles.linearSysNumOfVarEdit, 'string'))
dataInput = cell(rows, rows + 1)
%dataInput = cell(2,3)
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
%handles.matrixInput = get(tableInput,'data')
% hObject    handle to linearSysMatrixButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in linearSysSolveButton.
function linearSysSolveButton_Callback(hObject, eventdata, handles)
global tableInput
global rows
retrievedData = get(tableInput, 'data')
%= cell2mat(retrievedData)
augA = str2double(retrievedData)

sol = Gauss_Elim(augA)
% hObject    handle to linearSysSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Naive Gaussian Elimination
% x,y,z variables specification
% # of rows
% # of cols
% row input csv values?
function solutions = Gauss_Elim(augA)
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
n = length(augA) - 1
row = length(augA) - 1
col = row + 1

A = augA
A(:,col) = []
norminf = norm(A, inf);
norminf_inv = norm(inv(A), inf);

cond_num = norminf * norminf_inv;

iso_exp = floor(log10(cond_num*10));
fprintf('Error Magnification factors of the magnitude %d are possible.\n', iso_exp);
fprintf('Since Matlab defaults to double precision this means that \n');
fprintf('16 - %d = %d correct digits in the solution.\n', iso_exp, 16-iso_exp);

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
    solutions(q) = augA(q, row + 1)
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
toc;
fprintf('\n');

fprintf('Number of Operations = %d\n\n', opCount);

% -------------------------
% --- NONLINEAR SYSTEMS ---
% -------------------------


% --- Executes on button press in nonLinearSysSolveButton.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearSysSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonLinearSysSolveButton


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
