clear;
%% load data
data = load("datacell.mat");  % X: M*13 matrix. i step of void ratio(1), i step of stress(3), i step of strain(3),
data = data.datacell;        % i step of strain increment(3), i step of plastic strain(3)
datanum = 1:20;
M = (1:length(datanum))';     % Y: M*7 matrix. i+1 step of void ratio(1), i+1 step of stress(3), i+1 step of plastic starin increment(3)

%% set training, validation, and testing dataset
ratio = [0.8, 0.2];             %% ratio of training, validation, and testing dataset

[X, Y] = SetXY(data, datanum);
validation = DataValidation(X);
XData = (cell2mat(X))';
YData = (cell2mat(Y))';

%% Normalization
X_max=max(XData,[],2);X_min=min(XData,[],2);
Y_max=max(YData,[],2);Y_min=min(YData,[],2);

x = (XData-X_min)./(X_max-X_min);
y = (YData-Y_min)./(Y_max-Y_min);
%% divide training and testing sets
num_train = ceil(ratio(1)*length(y));
num_test = length(y) - num_train;
idx_train = randperm(length(y), num_train)';
idx_test = setdiff(1:length(y), idx_train);
    
xTrain = x(:,idx_train);XTrain = XData(:,idx_train);
yTrain = y(:,idx_train);YTrain = YData(:,idx_train);
xTest  = x(:,idx_test); XTest  = XData(:,idx_test);
yTest  = y(:,idx_test); YTest  = YData(:,idx_test);

f_input = size(xTrain, 1);
f_output = size(yTrain, 1);

ds = arrayDatastore([xTrain;yTrain]);
%% structure of LSTM
layers = [featureInputLayer(f_input)  
    fullyConnectedLayer(80)
    reluLayer          
    fullyConnectedLayer(80) 
    reluLayer     
    fullyConnectedLayer(80) 
    reluLayer 
    fullyConnectedLayer(f_output)];                
net = dlnetwork(layers);
net = dlupdate(@double,net);

%% hyperparameters 
InitialLearnRate = 0.001;
numEpochs = 5000;
miniBatchSize = 128;

decayRate = 0.001;

averageGrad = [];
averageSqGrad = [];

mbq = minibatchqueue(ds, ...
    MiniBatchSize=miniBatchSize, ...
    MiniBatchFormat="CB", ...
    OutputEnvironment="GPU");

accfun = dlaccelerate(@modelLoss);

figure(1)
clf
C = colororder;
Total_loss       = animatedline(Color=C(1,:),DisplayName="total loss");
e_loss        = animatedline(Color=C(2,:),DisplayName="e loss");
stress_loss        = animatedline(Color=C(4,:),DisplayName="stress loss");
pde_loss       = animatedline(Color=C(3,:),DisplayName="plastic strain loss");
legend
set(gca,"YScale",'log')
xlabel("epoch")
ylabel("Loss")
grid on
start=tic;

iteration = 0;
loss_record = [];
for epoch=1:numEpochs
    reset(mbq)

    while hasdata(mbq)
        iteration = iteration+1;
        xy = next(mbq);
        xt = xy(1:f_input,:);
        yt = xy(f_input+1:end,:);

        [loss, gradients,Y_pre,y_pre,y_pre_r] = ...
        dlfeval(dlaccelerate(@modelLoss), net, xt, yt, X_max,X_min,Y_max,Y_min);

        % learningRate = InitialLearnRate/(1+decayRate*epoch);
        learningRate = InitialLearnRate;
        
        [net, averageGrad, averageSqGrad] = adamupdate(net, gradients, averageGrad, averageSqGrad, iteration, learningRate);

        loss_record(end+1,:) = gather(extractdata(loss));
        
        addpoints(Total_loss, iteration, loss_record(iteration,1));
        addpoints(e_loss, iteration, loss_record(iteration,2));
        addpoints(stress_loss, iteration, loss_record(iteration,3));
        addpoints(pde_loss, iteration, loss_record(iteration,4));

        time = duration(0,0,toc(start),Format="hh:mm:ss");
        title("Iteration: " + iteration + ", Elapsed: " + string(time) + ", EE: GPU")
        drawnow
    end
end