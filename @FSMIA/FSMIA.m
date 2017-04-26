classdef FSMIA < handle
   properties
       filename
       Option
       Molecule
       Frame
       Result
       Intensity
   end
   
   methods
       function obj = FSMIA(filename)
           if nargin >0
               obj.filename = filename;
           end
           obj.Option = struct;
       end
       
       function obj = set.Option(obj,opt)
           obj.Option = opt;
       end
       
       function setoption(obj)
           opt = struct;
%            prompt = {'Threshold','Spot radius (pixels)','Pixel size (nm)',...
%                'Exclude region','Include only region',...
%                'Connect distance threshold (nm)','Fitting (fast or slow)',...
%                'Isolation Method (fast or slow)','Downsampling rate'...
%                'Illumination correction (on or off)','Background'...
%                'Wavelength (nm)','Numerical Aperture'};
            prompt = {'Spot radius (pixels)','Pixel size (nm)',...
               'Exclude region','Include only region',...
               'Connect distance threshold (nm)','Fitting (fast or slow)',...
               'Downsampling rate', 'Illumination correction (on or off)',...
               'Wavelength (nm)','Numerical Aperture','Threshold'};
           dlg_title = 'Set the options';
           def = {'5','160','0','0','320','slow','1','on','647','1.49','0.9'};
           answer = inputdlg(prompt,dlg_title,1,def);
           % opt.threshold = str2double(answer{1});
           opt.spotR = str2double(answer{1});
           opt.pixelSize = str2double(answer{2});
           exclude = str2double(answer{3});
           if isnan(exclude) || isempty(exclude)
               opt.exclude = false;
           else
               opt.exclude = exclude;
           end
           include = str2double(answer{4});
           if isnan(include) || isempty(include)
               opt.include = false;
           else
               opt.include = include;
           end
           opt.connectDistance = str2double(answer{5});
           opt.fitting = answer{6};
           % opt.isolation = answer{8};
           opt.ds = str2double(answer{7});
           opt.illumination = answer{8};
           % opt.bg = str2double(answer{11});
           opt.wavelength = str2double(answer{9});
           opt.na = str2double(answer{10});
           opt.threshold = str2double(answer{11});
           obj.Option = opt;
       end
   end
   
end
