clc
clear all
close all
%%

DirPath = uigetdir("");
cd(DirPath);

%%

data = load('dataset_yolov4');
data_table = data.data_table;
delete_index = [];
for i = 1 : height(data_table)
    a = data_table{i, 2};
    a = a{1};
    if isnan(a)
        delete_index = [delete_index i];
    end
end
data_table = removerows(data_table, 'ind', delete_index);

rng("default");
shuffled_indices = randperm(height(data_table));
indices = floor(0.7 * height(data_table));

training_data_indices = 1 : indices;
training_data_table = data_table(shuffled_indices(training_data_indices),:);

validation_data_indices = (indices + 1) : (indices + 1 + floor(0.2 * length(shuffled_indices)));
validation_data_table = data_table(shuffled_indices(validation_data_indices),:);

test_data_indices = validation_data_indices(end)+1 : length(shuffled_indices);
test_data_table = data_table(shuffled_indices(test_data_indices),:);

training_image_data_store = imageDatastore(training_data_table{:, 'image_locations'});
training_box_data_store = boxLabelDatastore(training_data_table(:, 'particle_positions'));

validation_image_data_store = imageDatastore(validation_data_table{:, 'image_locations'});
validation_box_data_store = boxLabelDatastore(validation_data_table(:, 'particle_positions'));

test_image_data_store = imageDatastore(test_data_table{:, 'image_locations'});
test_box_data_store = boxLabelDatastore(test_data_table(:, 'particle_positions'));

training_data = combine(training_image_data_store, training_box_data_store);
validation_data = combine(validation_image_data_store, validation_box_data_store);
test_data = combine(test_image_data_store, test_box_data_store);

input_size = [224 224 3];
number_classes = width(data_table) - 1;

rng("default")
trainingDataForEstimation = transform(training_data, @(data)preprocessData(data, input_size));
numAnchors = 9;
[anchors] = estimateAnchorBoxes(trainingDataForEstimation, numAnchors);

area = anchors(:, 1).*anchors(:,2);
[~,idx] = sort(area, "descend");

anchors = anchors(idx,:);
anchor_boxes = {anchors(1:3, :) ; anchors(4:6, :) ; anchors(7:9, :)};

class_name = "particle_positions";
detector = yolov4ObjectDetector("csp-darknet53-coco", class_name, anchor_boxes, InputSize = input_size);

augmentedTrainingData = transform(training_data, @augmentData);

options = trainingOptions("adam",...
    GradientDecayFactor=0.9,...
    SquaredGradientDecayFactor=0.999,...
    InitialLearnRate=0.001,...
    LearnRateSchedule="none",...
    MiniBatchSize=4,...
    L2Regularization=0.0005,...
    MaxEpochs=20,...
    BatchNormalizationStatistics="moving",...
    DispatchInBackground=true,...
    ResetInputNormalization=false,...
    Shuffle="every-epoch",...
    VerboseFrequency=20,...
    ValidationFrequency=100,...
    ValidationData=validation_data);

[detector, info] = trainYOLOv4ObjectDetector(augmentedTrainingData, detector, options);

preprocessed_test_data = transform(test_data, @(data)preprocessData(data, input_size));
detection_results = detect(detector, preprocessed_test_data);
[ap, recall, precision] = evaluateDetectionPrecision(detection_results, preprocessed_test_data);
                                                                                                                                                             
figure
plot(recall, precision)
xlabel('Recall')
ylabel('Precision')
grid on
title(sprintf('Average Precision = %.2f', ap))

detector_name = sprintf('detector-%s.mat', datetime('now', 'Format', 'yyyy-MM-dd HH;mm;ss'));
save(detector_name, 'detector');

function data = preprocessData(data, target_size)
    for ii = 1 : size(data, 1)
        I = data{ii, 1};
        image_size = size(I);
        
        bboxes = data{ii, 2};
    
        I = im2single(imresize(I, target_size(1:2)));
        scale = target_size(1:2) ./ image_size(1:2);
        bboxes = bboxresize(bboxes, scale);
        
        data(ii,1:2) = {I, bboxes};
    end
end

function data = augmentData(A)
    data = cell(size(A));
    for ii = 1:size(A,1)
        I = A{ii,1};
        bboxes = A{ii,2};
        labels = A{ii,3};
        sz = size(I);
    
        if numel(sz) == 3 && sz(3) == 3
            I = jitterColorHSV(I,...
                contrast=0.0,...
                Hue=0.1,...
                Saturation=0.2,...
                Brightness=0.2);
        end
        
        tform = randomAffine2d(XReflection=true,Scale=[1 1.1]);
        rout = affineOutputView(sz,tform,BoundsStyle="centerOutput");
        I = imwarp(I,tform,OutputView=rout);
        
        [bboxes,indices] = bboxwarp(bboxes,tform,rout,OverlapThreshold=0.25);
        labels = labels(indices);
        
        if isempty(indices)
            data(ii,:) = A(ii,:);
        else
            data(ii,:) = {I,bboxes,labels};
        end
    end
end