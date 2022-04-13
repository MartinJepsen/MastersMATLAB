clear all
damages = ones(6, 1);
damages(:, 2) = 0.85;
damages(:, 1) = 1:6;

for ijk = 1:size(damages, 1)
    dam = damages(ijk, :);
    estimate_models
    gaindesign_2
    gaindesign_3
    
end

%% Run single setup
dam = damages(1, :);
set_up
apply_gains