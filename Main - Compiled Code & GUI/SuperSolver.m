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
function SuperSolver_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SuperSolver (see VARARGIN)

% Choose default command line output for SuperSolver
handles.output = hObject;
%% Tabs Code
% Settings
TabFontSize = 10;
TabNames = {'Single Variable','Linear Systems','Nonlinear Systems'};
FigWidth = 0.5;

% Figure resize
set(handles.SimpleOptimizedTab,'Units','normalized')
pos = get(handles. SimpleOptimizedTab, 'Position');
set(handles. SimpleOptimizedTab, 'Position', [pos(1) pos(2) FigWidth pos(4)])

% Tabs Execution
handles = TabsFun(handles,TabFontSize,TabNames);

% Update handles structure
guidata(hObject, handles);
set(handles.singleVarOutputText, 'Max', 5);
set(handles.lSysOutputText, 'Max', 5);
set(handles.linearSysGuessButton, 'enable', 'off');
set(handles.omegaEdit, 'enable', 'off');
set(handles.singleVarx0Edit, 'enable', 'off');
set(handles.nonLinearSysBroydenMatrixEdit, 'enable', 'off')
warning('off', 'symbolic:sym:sym:DeprecateExpressions');
% UIWAIT makes SuperSolver wait for user response (see UIRESUME)
% uiwait(handles.SimpleOptimizedTab);

% --- TabsFun creates axes and text objects for tabs
function handles = TabsFun(handles,TabFontSize,TabNames)

% Set the colors indicating a selected/unselected tab
handles.selectedTabColor=get(handles.tab1Panel,'BackgroundColor');
handles.unselectedTabColor=handles.selectedTabColor-0.1;
handles.matrixFigure = {};
handles.guessFigure = {};
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

for i = 1:handles.TabsNumber
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
function varargout = SuperSolver_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%==============================
% Single Variable GUI elements
%==============================
% --- Executes on selection change in singleVarListBox.
function singleVarListBox_Callback(~, ~, handles)
% hObject    handle to singleVarListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
method = get(handles.singleVarListBox, 'value');
switch method
    case 1
        set(handles.bisectionAEdit, 'enable', 'on');
        set(handles.bisectionBEdit, 'enable', 'on');
        set(handles.singleVarx0Edit, 'enable', 'off');
    otherwise
        set(handles.bisectionAEdit, 'enable', 'off');
        set(handles.bisectionBEdit, 'enable', 'off');
        set(handles.singleVarx0Edit, 'enable', 'on');
end

% Hints: contents = cellstr(get(hObject,'String')) returns singleVarListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from singleVarListBox


% --- Executes during object creation, after setting all properties.
function singleVarListBox_CreateFcn(hObject, ~, ~)
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
function singleVarToleranceEdit_CreateFcn(hObject, ~, ~)
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
function singleVarIterationsEdit_CreateFcn(hObject, ~, ~)
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
function singleVarx0Edit_CreateFcn(hObject, ~, ~)
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
function singleVarFunctionEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to singleVarFunctionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function singleVarRootEdit_Callback(hObject, eventdata, handles)
% hObject    handle to singleVarRootEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of singleVarRootEdit as text
%        str2double(get(hObject,'String')) returns contents of singleVarRootEdit as a double


% --- Executes during object creation, after setting all properties.
function singleVarRootEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to singleVarRootEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bisectionAEdit_Callback(~, ~, ~)
% hObject    handle to bisectionAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bisectionAEdit as text
%        str2double(get(hObject,'String')) returns contents of bisectionAEdit as a double


% --- Executes during object creation, after setting all properties.
function bisectionAEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to bisectionAEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bisectionBEdit_Callback(~, ~, ~)
% hObject    handle to bisectionBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bisectionBEdit as text
%        str2double(get(hObject,'String')) returns contents of bisectionBEdit as a double


% --- Executes during object creation, after setting all properties.
function bisectionBEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to bisectionBEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in singleVarSolveButton.
function singleVarSolveButton_Callback(hObject, ~, handles)
% hObject    handle to singleVarSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Formatting text box
set(handles.singleVarOutputText, 'Max', 5);
set(handles.singleVarOutputText, 'string', '');
method = get(handles.singleVarListBox, 'value');

% Receiving input from GUI elements
syms infxn;
rawfxn = get(handles.singleVarFunctionEdit, 'string');
infxn = sym(rawfxn);
x0 = str2double(get(handles.singleVarx0Edit, 'string'));
Tol = str2double(get(handles.singleVarToleranceEdit, 'string'));
r = str2double(get(handles.singleVarRootEdit, 'string'));
maxIterations = str2double(get(handles.singleVarIterationsEdit, 'string'));



errFlag = false; % Indicates if the methods have stopped

% Checks for valid/complete input
if isempty(rawfxn) || isnan(x0) || isnan(Tol) || isnan(maxIterations) || abs(r) == Inf
    errorString = 'Invalid Input';
    set(handles.singleVarOutputText, 'string', errorString);
else
    % Empty values for later output
    set(handles.SimpleOptimizedTab, 'HandleVisibility', 'off');
    close all;
    set(handles.SimpleOptimizedTab, 'HandleVisibility', 'on');
    xList = zeros();
    errList = zeros();
    iterations = 0;
    opCount = 0;
    %timeStr;
    if isnan(r)
        r = Inf;    % Checks to see if user entered real rool, can run in different manner
    end
    outhandles = guidata(hObject);  % Passing handles to functions to work with textbox
    
    switch method
        case 1
            a = str2double(get(handles.bisectionAEdit, 'string'));
            b = str2double(get(handles.bisectionBEdit, 'string'));
            if isnan(a) || isnan(b)
                rangeErrorString='Missing range input';
                set(handles.singleVarOutputText, 'string', rangeErrorString);
                errFlag = true;
            elseif a > b
                rangeErrString2 = 'a must be less than b';
                set(handles.singleVarOutputText, 'string', rangeErrString2);
                errFlag = true;
            else
                [xList, errList, errFlag, iterations, opCount, timeStr] = Bisection(infxn, a, b, r, Tol, maxIterations, outhandles);
            end
        case 2
            [xList, errList, errFlag, iterations, opCount, timeStr] = Single_Var_FixedPoint(infxn, x0, r, Tol, maxIterations, outhandles);
        otherwise
            [xList, errList, errFlag, iterations, opCount, timeStr] = Single_Var_Newtons(infxn, x0, r, Tol, maxIterations, outhandles);
    end
    if errFlag == false
        calcError(infxn, matlabFunction(diff(infxn)), r, xList(end), outhandles, iterations, opCount, timeStr);
    end
    if length(xList) > 1
        f = figure(1);
        t = uitable(f);
        %figure
        iList = 0:iterations;
        subplot(1,2,1)
        pos = get(subplot(1,2,2),'position');
        delete(subplot(1,2,2));
        set(t, 'units', 'normalized');
        set(t, 'position', pos)
        set(t, 'ColumnWidth', {110 110});
        plot(iList, xList);
        xlabel('Iteration');
        ylabel('x');
        set(gca, 'xtick', 0:iterations+1);
        title('X Through Iterations');
        newRowIndices = strings;
        for j=0:iterations-1
            newRowIndices(j+1) = num2str(j);
        end
        newRowIndices(j+2) = num2str(j+1);
        set(t, 'RowName', newRowIndices);
        xStrings = cell(size(xList));
        for s = 1:iterations+1
            xStrings(s) = cellstr(num2str(xList(s), '%20.16f'));
        end
        if r ~= Inf
            errStrings = cell(size(errList));
            for s = 1:iterations+1
                errStrings(s) = cellstr(num2str(errList(s), '%20.16f'));
            end
            combinedList = [xStrings; errStrings]';
            set(t,'Data',combinedList); % Use the set command to change the uitable properties.
            set(t,'ColumnName',{'xi', 'Error'})
        else
            combinedList = xStrings';
            set(t,'Data',combinedList);
            set(t,'ColumnName',{'xi'})
            set(t, 'RowName', newRowIndices);
        end
    end
end

%====================================
% Linear Systems GUI Elements
%====================================

% --- Executes on selection change in linearSysListBox.
function linearSysListBox_Callback(~, ~, handles)
method = get(handles.linearSysListBox, 'value');
switch method
    case 1
        set(handles.linearSysGuessButton, 'enable', 'off');
        set(handles.omegaEdit, 'enable', 'off');
    case 2
        set(handles.linearSysGuessButton, 'enable', 'off');
        set(handles.omegaEdit, 'enable', 'off');
    case 3
        set(handles.linearSysGuessButton, 'enable', 'on');
        set(handles.omegaEdit, 'enable', 'off');
    otherwise
        set(handles.linearSysGuessButton, 'enable', 'on');
        set(handles.omegaEdit, 'enable', 'on');
end

% --- Executes during object creation, after setting all properties.
function linearSysListBox_CreateFcn(hObject, ~, ~)
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

function linearSysNumOfVarEdit_CreateFcn(hObject, ~, ~)
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
function linearSysIterEdit_CreateFcn(hObject, ~, ~)
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
function linearSysTolEdit_CreateFcn(hObject, ~, ~)
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
function omegaEdit_CreateFcn(hObject, ~, ~)
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
    columnNames = strings;
    for i = 1:rows
        columnNames(i) = ['x',num2str(i)];
    end
    columnNames(i+1) = 'b';
    figure,
    pos = get(gcf,'position');
    %set(gcf,'position',[pos(1:2) [660 120]])
    %Input table's creation
    tableInput = uitable('ColumnWidth',{50},...
        'Position',[0 0 pos(3) pos(4)],...% 600 80], ...
        'data',dataInput, ...
        'columnName', columnNames,...
        'ColumnEditable',true);
end

% --- Executes on button press in linearSysSolveButton.
function linearSysSolveButton_Callback(hObject, ~, handles)
linHandles = guidata(hObject);
emptyFlag = true;
global tableInput;
global initialMatrixInput;

try
    retrievedData = get(tableInput, 'data');
    emptyFlag = false;
catch
    set(handles.lSysOutputText, 'string', 'Error: Missing matrix input.');
    emptyFlag = true;
end

if emptyFlag == false
    augA = str2double(retrievedData);
    method = get(handles.linearSysListBox, 'value');
    Tol = str2double(get(handles.linearSysTolEdit, 'string'));
    omega = str2double(get(handles.omegaEdit, 'string'));
    iterations = str2double(get(handles.linearSysIterEdit, 'string'));
    set(handles.lSysOutputText, 'Max', 2);
    rows = str2double(get(handles.linearSysNumOfVarEdit, 'string'));
    
    if isnan(Tol) || isnan(iterations)
        fprintf('Need Input')
    else
        sol = [];
        errFlag = true;
        switch method
            case 1
                [sol, errFlag] = Gauss_Elim(augA, linHandles);
            case 2
                [sol, errFlag] = LU_Decomposition(augA, linHandles);
            case 3
                 try
                    retrievedGuess = get(initialMatrixInput, 'data');
                    errFlag = false;
                 catch
                    set(handles.lSysOutputText, 'string', 'Error: missing guess matrix input');
                     errFlag = true;
                 end
                if errFlag == false
                    P = str2double(retrievedGuess);
                    [sol, errFlag] = Jacobi(augA, P, Tol, iterations, linHandles);
                end
            otherwise
                if isnan(omega)
                    errStr = 'need input';
                    set(handles.lSysOutputText, 'string', errStr);
                else
                    try
                        retrievedGuess = get(initialMatrixInput, 'data');
                        errFlag = false;
                    catch
                        set(handles.lSysOutputText, 'string', 'Error: missing guess matrix input');
                        errFlag = true;
                    end
                    if errFlag == false
                        x0 = str2double(retrievedGuess);
                        [sol, errFlag] = SOR(augA, x0, omega, iterations, Tol, linHandles);
                    end
                end
        end
        set(handles.SimpleOptimizedTab, 'HandleVisibility', 'off');
        close all;
        set(handles.SimpleOptimizedTab, 'HandleVisibility', 'on');
        if errFlag ~= true
            rowNames = strings;
            for i = 1:rows
                rowNames(i) = ['x',num2str(i)];
            end
            f = figure;
            t = uitable(f);
            set(t, 'Data', sol');
            set(t, 'RowName', rowNames);
            set(t, 'ColumnName', 'Solutions');
        end
    end
end
% hObject    handle to linearSysSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in linearSysGuessButton.
function linearSysGuessButton_Callback(~, ~, handles)
set(handles.lSysOutputText, 'Max', 5);
global initialMatrixInput;

vars = str2double(get(handles.linearSysNumOfVarEdit, 'string'));

if isnan(vars) || vars < 1 || vars == Inf
    errorString = 'Invalid variable amount';
    set(handles.lSysOutputText, 'string', errorString);
else
    dataInput = cell(1, vars);
    figure,
    pos = get(gcf,'position');
    columnNames = strings;
    for i = 1:vars
        columnNames(i) = ['x',num2str(i)];
    end
    %set(gcf,'position',[pos(1:2) [660 120]])
    %Input table's creation
    initialMatrixInput = uitable('ColumnWidth',{70},...
        'Position',[0 0 pos(3) pos(4)],...
        'data',dataInput, ...
        'columnName', columnNames, ...
        'rowName', 'Guesses',...
        'ColumnEditable',true);
end

%==============================
% Linear Systems Methods
%==============================

% -------------------------
% --- NONLINEAR SYSTEMS ---
% -------------------------


% --- Executes on selection change in nonLinearSysListBox.
function nonLinearSysListBox_Callback(hObject, eventdata, handles)
method = get(handles.nonLinearSysListBox, 'value');
if method == 1
    set(handles.nonLinearSysBroydenMatrixEdit, 'enable', 'off')
else
    set(handles.nonLinearSysBroydenMatrixEdit, 'enable', 'on')
end
% hObject    handle to nonLinearSysListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nonLinearSysListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nonLinearSysListBox


% --- Executes during object creation, after setting all properties.
function nonLinearSysListBox_CreateFcn(hObject, eventdata, handles)
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nonLinearMultiVarSolveButton.
function nonLinearMultiVarSolveButton_Callback(hObject, eventdata, handles)
% hObject    handle to nonLinearMultiVarSolveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonLinearMultiVarSolveButton
method = get(handles.nonLinearSysListBox, 'value');
%disp(method);
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
% celldisp(ca);

%convert each char to sym type
vars = cell2sym(ca);
for i=1:length(variables)
    syms(variables(i))
end

%GET initial guess
x = str2num(get(handles.nonLinearSysInitialGuessEdit, 'string'));
 

%GET initial matrix
A = str2num(get(handles.nonLinearSysBroydenMatrixEdit, 'string'));

%GET system of equations
eqn_str = get(handles.nonLinearSysEqnListEdit, 'string');
eqns = sym(eqn_str);
% disp(eqns);

%GET actual root
root = str2num(get(handles.nonLinearSysRootEdit, 'string'));


%GET number of iterations
number_of_iterations = str2double(get(handles.nonLinearSySIterationsEdit, 'string'));
% disp(number_of_iterations);
nLSHandles = guidata(hObject);
switch method
    case 1
        [appr_x] = Multi_Var_Newton_Method(vars, x, eqns, number_of_iterations, nLSHandles);
        [forward, backward]=find_error(transpose(appr_x), root, eqns,vars);
        forwardSpec = sprintf('Forward error is %4.2f%s\n%s', forward);
        backwardSpec = sprintf('Backward error is %4.2f%s\n%s', backward);
       
        appr_x = sprintf('%0.5f, ',appr_x);
        appr_x = appr_x(1:end-2);
        
        answer = ['Approximate root: ', appr_x];
        string = sprintf('%s\n%s\n%s', answer, forwardSpec, backwardSpec);
        currStr = get(handles.nonLinearSysOutput, 'string');
        string = combineString(currStr, string);
        set(handles.nonLinearSysOutput, 'Max', 20, 'HorizontalAlignment', 'left', 'String',  string);
        
    case 2
        [appr_x] = Broyden_1_Method(x, vars, eqns, A, number_of_iterations);
        
        [forward, backward]=find_error(appr_x, root, eqns,vars);
        forwardSpec = sprintf('Forward error is %4.2f%s\n%s', forward);
        backwardSpec = sprintf('Backward error is %4.2f%s\n%s', backward);
        
        appr_x = sprintf('%0.5f, ',appr_x);
        appr_x = appr_x(1:end-2);
        
        answer = ['Approximate root: ', appr_x];
        string = sprintf('%s\n%s\n%s', answer, forwardSpec, backwardSpec);
        
        set(handles.nonLinearSysOutput, 'Max', 20, 'HorizontalAlignment', 'left', 'String',  string);
        
    otherwise
        [appr_x] = Broyden_2_Method(x, vars, eqns, A, number_of_iterations);
        [forward, backward]=find_error(appr_x, root, eqns,vars);
        forwardSpec = sprintf('Forward error is %4.2f%s\n%s', forward);
        backwardSpec = sprintf('Backward error is %4.2f%s\n%s', backward);
        
        appr_x = sprintf('%0.5f, ',appr_x);
        appr_x = appr_x(1:end-2);
        
        answer = ['Approximate root: ', appr_x];
        string = sprintf('%s\n%s\n%s', answer, forwardSpec, backwardSpec);
        set(handles.nonLinearSysOutput, 'Max', 20, 'HorizontalAlignment', 'left', 'String',  string);
end

% --- Executes during object creation, after setting all properties.
function nonLinearSysBroydenMatrixEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function lSysOutputText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lSysOutputText (see GCBO)

function nonLinearSysOutput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function nonLinearSysRootEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%========================
% Menu Bar Functionality
%========================

% --------------------------------------------------------------------
function FileMenuButton_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenuButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HelpMenuTag_Callback(hObject, eventdata, handles)
% hObject    handle to HelpMenuTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function exitMenuChoice_Callback(hObject, eventdata, handles)
close all;
% hObject    handle to exitMenuChoice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function svHelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to svHelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpString = sprintf('Input values into their corresponding fields.\nIf the real root is unknown, it can be left blank.\nThen click solve and the selected algorithm will run.');
helpdlg(helpString,...
    'Single Variables');


% --------------------------------------------------------------------
function lsHelpMenu_Callback(~, ~, ~)
% hObject    handle to lsHelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpString = sprintf('Input values into their corresponding fields.\nThen generate a matrix to fill with the coefficients of the function based on the number of variables specified.\nFor Jacobi and SOR, generate and fill the inital vector in the same manner.');
helpdlg(helpString,...
    'Linear Systems');


% --------------------------------------------------------------------
function nlsHelpMenu_Callback(hObject, eventdata, handles)
% hObject    handle to nlsHelpMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpString = sprintf('Input values into their corresponding fields.\nInput a list of variable identifers that will be used in the functions.\nAfterwards, input the functions and other values using Matlab Matrix notation.\nInput must be in square brackets []\nInput is done in rows with elements separated by a space.\nRows are then separated by a semicolon.');
helpdlg(helpString,...
    'Noninear Systems');

% --------------------------------------------------------------------
function aboutMenu_Callback(hObject, eventdata, handles)
% hObject    handle to aboutMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aboutString = sprintf('This program calculates roots and solutions to functions using the following methods:\nSingle Variables: Bisection, Fixed-Point Iteration, Newton''s\nLinear Systems:Gaussian Elimination, LU Decomposition, Jacobi, Successive-Over-Relaxation\nNonlinear Systems: Multivariate Newton''s, Broyden I, BroydenII.\n\nContributors to this program are:\nApril Hudspeth\nSarah Jansen\nLatisha Konz\nLong Truong\nRui Zhang');
helpdlg(aboutString,...
    'Noninear Systems');
