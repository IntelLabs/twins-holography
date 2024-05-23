function [I] = getimage(Peg)
    % Peg in camera parameters [Exposure Gain] OR
    % cell {[Exposure Gain] vid src}
    if(iscell(Peg))
        eg = Peg{1};
        vid = Peg{2};
        src = Peg{3};
        
        if numel(eg)>0
            src.ExposureAuto = 'Off';
            src.ExposureTime = eg(1);
        end
        if numel(eg)>1
            src.GainAuto = 'Off';
            src.Gain = eg(2);
        end
        if numel(eg)>2
            src.AcquisitionFrameRate = eg(3);
        end

        I = getsnapshot(vid); 
    else
        vid = cam_init(Peg);

        I = getsnapshot(vid);        

        delete(vid);
    end
end

