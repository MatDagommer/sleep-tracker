classdef ViewerDreemRecord < handle
    
    properties(Access=protected)
        
    end

    properties
        after_disp          % time period to display after marker 
        ax_list             % axes of the figure
        before_disp         % time period to display before marker
        color_list          % list of colors used in the figure
        duration            % duration of the signals = number of point
        fig                 % figure
        fs                  % sampling frequency
        ind                 % index of the current time marker, used to choose which part of signal to plot
        nb_channel          % number of channels
        nb_ind              % number of time markers
        nextbutton          % button '>'
        prevbutton          % button '<'
        rec_colors          % colors of the rectangles
        rec_durations       % duration of the rectangles  
        rec_starts          % begining of the rectangles
        signals             % signals to plot
        plot_colors         % colors of the signals
        slidind             % slider
        time_events         % markers of events
        time_views          % markers of visualisation
        titles              % graph titles
        verticals           % timestamps of vertical lines
        vertical_colors     % colors of vertical lines
        ymins               % ymin for each graph 
        ymaxs               % ymax for each plot
    end
    
    %methods
    methods
        
        %% Constructor
        function obj = ViewerDreemRecord(signals, varargin)
            % Check number of parameters
            if nargin < 1 || mod(length(varargin),2) ~= 0
              error('Incorrect number of parameters.');
            end
            
            % Parse parameter list
            for i = 1:2:length(varargin)
                if ~ischar(varargin{i})
                    error(['Parameter ' num2str(i+2) ' is not a property.']);
                end
                switch(lower(varargin{i}))
                    case 'fs'
                        obj.fs = varargin{i+1};
                        if ~isint(obj.fs) || obj.fs <=0
                            error('Incorrect value for property ''Fs''.');
                        end
                    case 'displaywindow'
                        window_display = varargin{i+1};
                        if ~isfloat(window_display) || length(window_display)~=2
                            error('Incorrect value for property ''DisplayWindow''.');
                        end
                    otherwise
                        error(['Unknown property ''' num2str(varargin{i}) '''.']);
                end
            end
            
            %init variables
            for ch=1:length(signals)
                if size(signals{ch},1)>size(signals{ch},2)
                    signals{ch} = signals{ch}';
                end
            end
            obj.signals = signals;
            
            if isempty(obj.fs)
                obj.fs = 250;
            end
            if exist('window_display', 'var')
                obj.before_disp = window_display(1);
                obj.after_disp = window_display(2);
            else
                obj.before_disp = 3;
                obj.after_disp = 4;
            end
            obj.nb_channel = max(size(obj.signals));
            obj.duration = max(cellfun(@(v)max(v(1,:)), obj.signals));
            obj.titles =  cell(1,obj.nb_channel);
            obj.ymins = -5000 * ones(1,obj.nb_channel);
            obj.ymaxs = 5000 * ones(1,obj.nb_channel);
            obj.ind = 1;
            obj.time_views = 0:1:obj.duration;
            obj.time_events = [];
            obj.nb_ind = length(obj.time_views);
            
            obj.verticals = cell(1,obj.nb_channel);
            obj.vertical_colors = cell(1,obj.nb_channel);
            obj.rec_starts = cell(1,obj.nb_channel);
            obj.rec_durations = cell(1,obj.nb_channel);
            obj.rec_colors = cell(1, obj.nb_channel);
            obj.color_list = cell(0);
            obj.plot_colors = {'y','k','b','r','g',[0.5 0.5 0.5]};
        end
        
        %set titles
        function set_title(obj, title, channel)
            obj.titles{channel} = title;
        end
        
        % add verticals bar
        function add_verticals(obj, timestamps, channel, color)
            if nargin < 4
                color = 'b';
            end
            ind_color = length(obj.color_list) + 1;
            obj.color_list{ind_color} = color;
            obj.verticals{channel} = [obj.verticals{channel} timestamps];
            obj.vertical_colors{channel} = [obj.vertical_colors{channel} ones(1,length(timestamps))*ind_color];
        end
        
        % add rectangles
        function add_rectangles(obj, timestamps, durations, channel, color)
            if nargin < 5
                color = [0.75 0.75 0.75];
            end
            if length(timestamps)~=length(durations)
                error('start and durations must have th same length')
            end
            ind_color = length(obj.color_list) + 1;
            obj.color_list{ind_color} = color;
            obj.rec_starts{channel}   = [obj.rec_starts{channel} timestamps];
            obj.rec_durations{channel} = [obj.rec_durations{channel} durations];
            obj.rec_colors{channel}    = [obj.rec_colors{channel} ones(1,length(timestamps))*ind_color];
        end
        
        % change time views
        function set_index_views(obj, timestamps)
            obj.time_views = sort(timestamps);
            obj.ind = 1;
            obj.nb_ind = length(timestamps);
        end
        function reset_time_views(obj)
            obj.time_views = 0:1:obj.duration;
            obj.ind = 1;
            obj.nb_ind = length(obj.time_views);
        end
        % change time events
        function set_time_events(obj, timestamps)
            obj.time_events = sort(timestamps);
        end
        
        % change ylimits
        function set_ymin(obj, channel, ymin)
            obj.ymins(channel) = ymin;
        end
        function set_ymax(obj, channel, ymax)
            obj.ymaxs(channel) = ymax;
        end
        %autoscale ylimits
        function autoscale_ylim(obj, channel)
            std_sig = std(obj.signals{channel}(2,:));
            obj.ymins(channel) = - std_sig * 5;
            obj.ymaxs(channel) = std_sig * 5;
        end
        
        
        %% Create the figure
        function run_window(obj)
            
            %  Create and then hide the UI as it is being constructed.
            obj.fig = figure('Name','SignalObserver','Visible','off','units','normalized','outerposition',[0 0 1 1]);

            % Construct the components.
            obj.prevbutton = uicontrol('Style','pushbutton', 'String','<','units','normalized','Position',[0.80 0.03 0.05 0.045], 'Callback',@(handle,event)prevbutton_Callback(obj,handle,event));
            obj.nextbutton = uicontrol('Style','pushbutton', 'String','>','units','normalized','Position',[0.85 0.03 0.05 0.045], 'Callback',@(handle,event)nextbutton_Callback(obj,handle,event));
            obj.slidind = uicontrol('Style','slider','units','normalized', 'Position',[0.08 0.03 0.55 0.025], 'value', obj.ind, 'min',1, 'max',obj.nb_ind,...
                'Callback', @(handle,event)slidind_Callback(obj,handle,event));
            
            set(obj.fig,'KeyPressFcn', @(handle,event)key_pressed_fcn(obj,handle,event));
            
            obj.ax_list = cell(1,obj.nb_channel);
            for ch=1:obj.nb_channel
                obj.ax_list{ch} = subplot(obj.nb_channel,1,ch);
            end
            
            center = obj.time_views(obj.ind);
            if center > obj.before_disp
                x0 = center - obj.before_disp;
            else
                x0 = 0;
            end
            if center + obj.after_disp <= obj.duration
                xf = center + obj.after_disp;
            else
                xf = obj.duration;
            end
            
            %Plot
            for ch=1:obj.nb_channel
                set(obj.fig,'CurrentAxes',obj.ax_list{ch});
                xdata = obj.signals{ch}(1,:);
                t = xdata(xdata>=x0 & xdata<=xf);
                baseline = zeros(1, length(t));
                
                for k=2:size(obj.signals{ch},1)
                    ydata = obj.signals{ch}(k,xdata>=x0 & xdata<=xf);
                    plot(t,ydata,'LineWidth',1,'Color',obj.plot_colors{k}), hold on,
                end
                
                plot(t,baseline,'LineWidth',0.3, 'Color', [0.7 0.7,0.7]),
                y_min = obj.ymins(ch);
                y_max = obj.ymaxs(ch);
                ylim([y_min y_max]);
                xlim([x0 xf]);
                title(obj.ax_list{ch}, obj.titles{ch});
                %add verticals to plot
                if ~isempty(obj.verticals{ch})
                    vert = obj.verticals{ch};
                    vert_color = obj.vertical_colors{ch};
                    
                    vert_ok = vert>=x0 & vert<=xf;
                    vert = vert(vert_ok);
                    vert_color = vert_color(vert_ok);
                    
                    for i=1:length(vert)
                        hold on, line([vert(i) vert(i)], [y_min y_max], 'Color', obj.color_list{vert_color(i)});
                    end
                end
                
                %add rectangles to plot
                if ~isempty(obj.rec_starts{ch})
                    rec_st = obj.rec_starts{ch};
                    rec_dur = obj.rec_durations{ch};
                    rec_color = obj.rec_colors{ch};
                    
                    rec_ok = rec_st>=x0 & rec_st<xf;
                    rec_st = rec_st(rec_ok);
                    rec_dur = rec_dur(rec_ok);
                    rec_color = rec_color(rec_ok);
                    
                    for i=1:length(rec_st)
                        if rec_st(i) + rec_dur(i) > xf
                            rec_dur(i) = xf - rec_st(i);
                        end
                        x_rec = [rec_st(i), rec_st(i)+rec_dur(i), rec_st(i)+rec_dur(i), rec_st(i)];
                        y_rec = [y_min, y_min, y_max, y_max];
                        hold on, patch(x_rec, y_rec, obj.color_list{rec_color(i)}, 'EdgeColor', 'None', 'FaceAlpha', 0.2); 
                    end
                end
            end

            % Make the window visible.
            obj.fig.Visible = 'on';    
        end
        
        
        %% Events and Callback
        function key_pressed_fcn(obj,handle,event)
            switch get(obj.fig, 'CurrentKey')
                case 'leftarrow'
                    prev(obj)
                case 'rightarrow'
                    next(obj)
                case 'pagedown'
                    prev(obj,10)
                case 'pageup'
                    next(obj,10)
                case 'downarrow'
                    prev_event(obj)
                case 'uparrow'
                    next_event(obj)
            end
        end
        
        function slidind_Callback(obj,handle,event)
            obj.ind = floor(handle.Value);
            if obj.ind==0
               obj.ind = obj.nb_ind;
            end
            plot_curve(obj)
        end
        
        function prevbutton_Callback(obj,handle,event)
           prev(obj); 
        end
        
        function prev(obj, steps)
           if nargin < 2
               steps = 1;
           end
           obj.ind = mod((obj.ind - steps), obj.nb_ind);
           if obj.ind==0
               obj.ind = obj.nb_ind;
           end
           plot_curve(obj)
        end
        
        function nextbutton_Callback(obj,handle,event)
           next(obj); 
        end
        
        function next(obj, steps)
           if nargin < 2
               steps = 1;
           end
           obj.ind = mod((obj.ind + steps), obj.nb_ind);
           if obj.ind==0
               obj.ind = obj.nb_ind;
           end
           plot_curve(obj)
        end
        
        function prev_event(obj)
            if ~isempty(obj.time_events)
                current_time = obj.time_views(obj.ind);
                if any(obj.time_events < current_time)
                    tmp_ev = obj.time_events(find(obj.time_events < current_time, 1, 'last'));
                    [~, idx] = min(abs(obj.time_views-tmp_ev));
                    obj.ind = idx;
                    plot_curve(obj)
                end
            end
        end
        function next_event(obj)
            if ~isempty(obj.time_events)
                current_time = obj.time_views(obj.ind);
                if any(obj.time_events > current_time)
                    tmp_ev = obj.time_events(find(obj.time_events > current_time, 1, 'first'));
                    [~, idx] = min(abs(obj.time_views-tmp_ev));
                    obj.ind = idx;
                    plot_curve(obj)
                end
            end
        end
        
        %% Update curves
        function plot_curve(obj)
            center = obj.time_views(obj.ind);
            if center > obj.before_disp + 1
                x0 = center - obj.before_disp;
            else
                x0 = 0;
            end
            if center + obj.after_disp <= obj.duration
                xf = center + obj.after_disp;
            else
                xf = obj.duration;
            end
            
            %Plot
            for ch=1:obj.nb_channel
                set(obj.fig,'CurrentAxes',obj.ax_list{ch});
                cla;
                xdata = obj.signals{ch}(1,:);
                t = obj.signals{ch}(1,xdata>=x0 & xdata<=xf);
                baseline = zeros(1, length(t)); 
                
                for k=2:size(obj.signals{ch},1)
                    ydata = obj.signals{ch}(k,xdata>=x0 & xdata<=xf);
                    plot(t,ydata,'LineWidth',1,'Color',obj.plot_colors{k}), hold on,
                end
                
                plot(t,baseline,'LineWidth',0.3, 'Color', [0.7 0.7,0.7]), 
                y_min = obj.ymins(ch);
                y_max = obj.ymaxs(ch);
                ylim([y_min y_max]);
                xlim([x0 xf]);
                
                %add verticals to plot
                if ~isempty(obj.verticals{ch})
                    vert = obj.verticals{ch};
                    vert_color = obj.vertical_colors{ch};
                    
                    vert_ok = vert>=x0 & vert<=xf;
                    vert = vert(vert_ok);
                    vert_color = vert_color(vert_ok);
                    
                    for i=1:length(vert)
                        hold on, line([vert(i) vert(i)], [y_min y_max], 'Color', obj.color_list{vert_color(i)});
                    end
                end
                
                %add rectangles to plot
                if ~isempty(obj.rec_starts{ch})
                    rec_st = obj.rec_starts{ch};
                    rec_dur = obj.rec_durations{ch};
                    rec_color = obj.rec_colors{ch};
                    
                    rec_ok = rec_st>=x0 & rec_st<xf;
                    rec_st = rec_st(rec_ok);
                    rec_dur = rec_dur(rec_ok);
                    rec_color = rec_color(rec_ok);
                    
                    for i=1:length(rec_st)
                        if rec_st(i) + rec_dur(i) > xf
                            rec_dur(i) = xf - rec_st(i);
                        end
                        x_rec = [rec_st(i), rec_st(i)+rec_dur(i), rec_st(i)+rec_dur(i), rec_st(i)];
                        y_rec = [y_min, y_min, y_max, y_max];
                        hold on, patch(x_rec, y_rec, obj.color_list{rec_color(i)}, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
                    end                    
                end
            end
            
        end
          
        
    end
        
             
end