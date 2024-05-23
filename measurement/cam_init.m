function [vid, src] = cam_init(Peg)
    vid = videoinput('gentl', 1, 'Mono12');
    vid.FramesPerTrigger = 1;
    if numel(Peg)>3
        triggerconfig(vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
    end
    src = getselectedsource(vid);
    if nargin > 0
        if numel(Peg)>0
            src.ExposureAuto = 'Off';
            src.ExposureTime = Peg(1);
        end
        if numel(Peg)>1
            src.GainAuto = 'Off';
            src.Gain = Peg(2);
        end
        if numel(Peg)>2
            src.AcquisitionFrameRate = Peg(3);
            src.AcquisitionFrameRateEnable = 'True';
        end
        if numel(Peg)>3 % synchronize camera
            src.LineSelector = 'Line2';
            src.LineMode = 'Input';
            src.TriggerSource = 'Line2';
            src.TriggerSelector = 'FrameStart';
            src.TriggerActivation = 'RisingEdge';
            src.TriggerDelay = Peg(4);
            src.TriggerMode = 'On';
        end
    end
end
