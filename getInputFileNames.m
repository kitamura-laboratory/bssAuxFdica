function [dirPath,fileName] = getInputFileNames(dataNo)

switch(dataNo)
    case 1
        dirPath = "./dataset/dev1_female4_src_12_E2A_conv/";
        fileName(1) = "dev1_female4_src_1_E2A_pos050130_mic2123_conv.wav";
        fileName(2) = "dev1_female4_src_2_E2A_pos050130_mic2123_conv.wav";
    case 2
        dirPath = "./dataset/dev1_female4_src_34_E2A_conv/";
        fileName(1) = "dev1_female4_src_3_E2A_pos050130_mic2123_conv.wav";
        fileName(2) = "dev1_female4_src_4_E2A_pos050130_mic2123_conv.wav";
    case 3
        dirPath = "./dataset/dev1_male4_src_12_E2A_conv/";
        fileName(1) = "dev1_male4_src_1_E2A_pos050130_mic2123_conv.wav";
        fileName(2) = "dev1_male4_src_2_E2A_pos050130_mic2123_conv.wav";
    case 4
        dirPath = "./dataset/dev1_male4_src_34_E2A_conv/";
        fileName(1) = "dev1_male4_src_3_E2A_pos050130_mic2123_conv.wav";
        fileName(2) = "dev1_male4_src_4_E2A_pos050130_mic2123_conv.wav";
    case 5
        dirPath = "./dataset/dev1_female4_src_12_JR2_conv/";
        fileName(1) = "dev1_female4_src_1_JR2_pos060120_mic2123_conv.wav";
        fileName(2) = "dev1_female4_src_2_JR2_pos060120_mic2123_conv.wav";
    case 6
        dirPath = "./dataset/dev1_female4_src_34_JR2_conv/";
        fileName(1) = "dev1_female4_src_3_JR2_pos060120_mic2123_conv.wav";
        fileName(2) = "dev1_female4_src_4_JR2_pos060120_mic2123_conv.wav";
    case 7
        dirPath = "./dataset/dev1_male4_src_12_JR2_conv/";
        fileName(1) = "dev1_male4_src_1_JR2_pos060120_mic2123_conv.wav";
        fileName(2) = "dev1_male4_src_2_JR2_pos060120_mic2123_conv.wav";
    case 8
        dirPath = "./dataset/dev1_male4_src_34_JR2_conv/";
        fileName(1) = "dev1_male4_src_3_JR2_pos060120_mic2123_conv.wav";
        fileName(2) = "dev1_male4_src_4_JR2_pos060120_mic2123_conv.wav";
end
end