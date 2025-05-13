clc
clear all
close all
%%

DirPath = uigetdir("");
cd(DirPath);

DirInfo=dir;
[NumOfFile,Dummy]=size(DirInfo);

FileList_png = [];
FileList_txt = [];

FileCount_png = 0;
FileCount_txt = 0;

for n = 1 : NumOfFile
    if(DirInfo(n).isdir == 0)
       FileName=DirInfo(n).name;
       Ext = FileName(end-3:end);
       
       if ((strcmp(Ext,'.png') == 1))
           FileCount_png = FileCount_png +1;
           FileList_png{FileCount_png,1} = FileName;
           disp(FileName);
       end
       
       if ((strcmp(Ext,'.txt') == 1) && contains(FileName,'info')==0)
           FileCount_txt = FileCount_txt +1;
           FileList_txt{FileCount_txt,1} = FileName;
           disp(FileName);
       end
       
    end
end

FileList_png = natsortfiles(FileList_png);
FileList_txt = natsortfiles(FileList_txt);



if(FileCount_png == FileCount_txt)
    data_size = FileCount_png;
else
    data_size = 0;
    fprintf('Number of file not same\n');
end

data_table = table('Size', [data_size 2], 'VariableTypes', {'cellstr', 'cell'}, 'VariableNames', {'image_locations', 'particle_positions'});
for n = 1:data_size
    FileName = sprintf('%s',FileList_png{n,1});
    data_table(n, 1) = cellstr(FileName);
end

for n = 1:data_size
    FileName_txt = sprintf('%s',FileList_txt{n,1});
    data = readmatrix(FileName_txt);
    data_table{n, 2} = {data};
end

save('dataset_yolov4.mat', 'data_table');


%%
FRAME_RATE = 10;
MovieFileName = sprintf('Learning data movie.avi');
myVideo = VideoWriter(MovieFileName);    
myVideo.FrameRate = FRAME_RATE;        
myVideo.Quality = 100;                 
open(myVideo);
MOVIE_BUFFER = [];

FOV_HEIGHT = 128;
FOV_WIDTH = 128;

for n = 1:data_size
    ImageFileName = sprintf('%s',FileList_png{n,1});
    Image = imread(ImageFileName);

    [Height,Width,Depth] = size(Image);
    
    FileName_txt = sprintf('%s',FileList_txt{n,1});
    
    data = [];
    data = readmatrix(FileName_txt);
    
    if(isnan(data) == 0)
        ImageBoxMark = insertShape(Image,'rectangle', data);
    else
        ImageBoxMark = Image;
    end
    
    Image = insertShape(Image,'Rectangle',[1  1  Width Height],'LineWidth', 1,'Color',[68 146 197]);
    ImageBoxMark = insertShape(ImageBoxMark,'Rectangle',[1  1  Width Height],'LineWidth', 1,'Color',[100 189 101]);
        
    if(Height == Width)
        MOVIE_BUFFER(1:Height, Width*0+1 :Width*1, : ) = (Image);  
        MOVIE_BUFFER(1:Height, Width*1+1 :Width*2, : ) = (ImageBoxMark);    
    elseif(Height ~= Width)
        MOVIE_BUFFER(Height*0+1 :Height*1, 1:Width,: ) = (Image);  
        MOVIE_BUFFER(Height*1+1 :Height*2, 1:Width,: ) = (ImageBoxMark);    
    end
    writeVideo(myVideo,uint8(MOVIE_BUFFER));
    
    YS_OUT_IMAGE = Height/2 - round(FOV_HEIGHT/2)+1;
    YE_OUT_IMAGE = Height/2 + round(FOV_HEIGHT/2);
    XS_OUT_IMAGE = Width/2 - round(FOV_WIDTH/2)+1;
    XE_OUT_IMAGE = Width/2 + round(FOV_WIDTH/2);
    CROP_OUT_IMAGE = ImageBoxMark(YS_OUT_IMAGE: YE_OUT_IMAGE, XS_OUT_IMAGE: XE_OUT_IMAGE,:);

end

close(myVideo);




