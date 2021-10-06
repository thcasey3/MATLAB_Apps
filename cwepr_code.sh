#!/bin/sh

#  cwepr_code.sh
#  
#
#  Created by Thomas Casey on 10/6/21.
#  
classdef cwEPR < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        cwEPRapp                     matlab.ui.Figure
        UIAxes                       matlab.ui.control.UIAxes
        OpenDataButton               matlab.ui.control.Button
        SysTabGroup                  matlab.ui.container.TabGroup
        Sys1Tab                      matlab.ui.container.Tab
        Sys1Switch                   matlab.ui.control.Switch
        Sys1Edit                     matlab.ui.control.EditField
        Sys1ListBox                  matlab.ui.control.ListBox
        Sys2Tab                      matlab.ui.container.Tab
        Sys2Switch                   matlab.ui.control.Switch
        Sys2Edit                     matlab.ui.control.EditField
        Sys2ListBox                  matlab.ui.control.ListBox
        Sys3Tab                      matlab.ui.container.Tab
        Sys3Switch                   matlab.ui.control.Switch
        Sys3Edit                     matlab.ui.control.EditField
        Sys3ListBox                  matlab.ui.control.ListBox
        Sys4Tab                      matlab.ui.container.Tab
        Sys4Switch                   matlab.ui.control.Switch
        Sys4Edit                     matlab.ui.control.EditField
        Sys4ListBox                  matlab.ui.control.ListBox
        ExpTab                       matlab.ui.container.Tab
        ExpEdit                      matlab.ui.control.EditField
        FromDataButton               matlab.ui.control.Button
        ExpListBox                   matlab.ui.control.ListBox
        FrequencySliderLabel         matlab.ui.control.Label
        FrequencySlider              matlab.ui.control.Slider
        StaticExpListBox             matlab.ui.control.ListBox
        FreqSpinner                  matlab.ui.control.Spinner
        OptTab                       matlab.ui.container.Tab
        OptEdit                      matlab.ui.control.EditField
        OptListBox                   matlab.ui.control.ListBox
        ProcessingTabGroup           matlab.ui.container.TabGroup
        BaselineTab                  matlab.ui.container.Tab
        LoadButton                   matlab.ui.control.Button
        BackgroundFilenameLabel      matlab.ui.control.Label
        SpectrumUsageSliderLabel     matlab.ui.control.Label
        SpectrumUsageSlider          matlab.ui.control.Slider
        PolynomialFitButton          matlab.ui.control.Button
        ApplyButton                  matlab.ui.control.Button
        BackgroundButton             matlab.ui.control.Button
        RevertBaseButton             matlab.ui.control.Button
        fromendsLabel                matlab.ui.control.Label
        msbackadjButton              matlab.ui.control.Button
        RequiresBioinformaticsToolboxLabel  matlab.ui.control.Label
        NormalizeBaseCheckBox        matlab.ui.control.CheckBox
        OrderSpinnerLabel            matlab.ui.control.Label
        PolyOrderSpinner             matlab.ui.control.Spinner
        SmoothTab                    matlab.ui.container.Tab
        SmoothingOptionsButtonGroup  matlab.ui.container.ButtonGroup
        MovingAverageButton          matlab.ui.control.RadioButton
        BinomialButton               matlab.ui.control.RadioButton
        FlatButton                   matlab.ui.control.RadioButton
        SavitskyGolayButton          matlab.ui.control.RadioButton
        wdenoiseButton               matlab.ui.control.RadioButton
        pdifEditFieldLabel           matlab.ui.control.Label
        pdifEditField                matlab.ui.control.EditField
        LevelLabel                   matlab.ui.control.Label
        wdOptionsEditField           matlab.ui.control.EditField
        RequiresWavletToolboxLabel   matlab.ui.control.Label
        orfloorlog2nPointsLabel      matlab.ui.control.Label
        SmoothingBusyLabel           matlab.ui.control.Label
        StepsSpinnerLabel            matlab.ui.control.Label
        StepsSpinner                 matlab.ui.control.Spinner
        ShowButton                   matlab.ui.control.Button
        SmoothingApplyButton         matlab.ui.control.Button
        SmoothingRevertButton        matlab.ui.control.Button
        IntegrateTab                 matlab.ui.container.Tab
        IntegrateButton              matlab.ui.control.Button
        ShowIntButtonGroup           matlab.ui.container.ButtonGroup
        SingleButton                 matlab.ui.control.RadioButton
        DoubleButton                 matlab.ui.control.RadioButton
        BothButton                   matlab.ui.control.RadioButton
        NeitherButton                matlab.ui.control.RadioButton
        BaselineCorrButtonGroup      matlab.ui.container.ButtonGroup
        NoneButton                   matlab.ui.control.RadioButton
        PolynomialorderButton        matlab.ui.control.RadioButton
        IntOrderSpinner              matlab.ui.control.Spinner
        IntegralTextAreaLabel        matlab.ui.control.Label
        DIEditField                  matlab.ui.control.TextArea
        DoubleLabel                  matlab.ui.control.Label
        ndMomentButton               matlab.ui.control.Button
        TrapzCheckBox                matlab.ui.control.CheckBox
        SumCheckBox                  matlab.ui.control.CheckBox
        SecondMomentEdit             matlab.ui.control.EditField
        IntSpectrumUsageSlider       matlab.ui.control.Slider
        SpectrumUsagefromendsLabel   matlab.ui.control.Label
        AutoButton                   matlab.ui.control.Button
        AutoBusyLabel                matlab.ui.control.Label
        QuantifyTab                  matlab.ui.container.Tab
        Step1Label                   matlab.ui.control.Label
        Step2Label                   matlab.ui.control.Label
        DIEditFieldLabel             matlab.ui.control.Label
        aEditField                   matlab.ui.control.NumericEditField
        x2Label                      matlab.ui.control.Label
        Label                        matlab.ui.control.Label
        bEditField                   matlab.ui.control.NumericEditField
        xLabel                       matlab.ui.control.Label
        EditField_2Label             matlab.ui.control.Label
        cEditField                   matlab.ui.control.NumericEditField
        Step2Label_2                 matlab.ui.control.Label
        Step4Label                   matlab.ui.control.Label
        FindconcentrationButton      matlab.ui.control.Button
        Step3Label                   matlab.ui.control.Label
        SpinCEditField               matlab.ui.control.EditField
        DIfactorEditField            matlab.ui.control.EditField
        hintsLabel                   matlab.ui.control.Label
        hint2Label                   matlab.ui.control.Label
        AdjustTab                    matlab.ui.container.Tab
        AddNoiseButton               matlab.ui.control.Button
        RevertAdjustmentsButton      matlab.ui.control.Button
        NoiseModelButtonGroup        matlab.ui.container.ButtonGroup
        fButton                      matlab.ui.control.RadioButton
        UniformButton                matlab.ui.control.RadioButton
        GaussianButton               matlab.ui.control.RadioButton
        ApplyAdjustmentsButton       matlab.ui.control.Button
        InterpolateButton            matlab.ui.control.Button
        NoisySimButton               matlab.ui.control.Button
        NoiseDTADSCCheckBox          matlab.ui.control.CheckBox
        NoiseasciiCheckBox           matlab.ui.control.CheckBox
        NoisetoworkspaceCheckBox     matlab.ui.control.CheckBox
        NoisetodataCheckBox          matlab.ui.control.CheckBox
        SNRSpinnerLabel              matlab.ui.control.Label
        SNRSpinner                   matlab.ui.control.Spinner
        PointsSpinnerLabel           matlab.ui.control.Label
        PointsSpinner                matlab.ui.control.Spinner
        ToolsTab                     matlab.ui.container.Tab
        ExportDataButton             matlab.ui.control.Button
        ConstantsDropDownLabel       matlab.ui.control.Label
        ConstantsDropDown            matlab.ui.control.DropDown
        ConversionToolButton         matlab.ui.control.Button
        PeriodicTableButton          matlab.ui.control.Button
        ExportDTADSCCheckBox         matlab.ui.control.CheckBox
        ExportasciiCheckBox          matlab.ui.control.CheckBox
        ExporttoworkspaceCheckBox    matlab.ui.control.CheckBox
        ShowButtonGroup              matlab.ui.container.ButtonGroup
        DataOnlyButton               matlab.ui.control.RadioButton
        SimOnlyButton                matlab.ui.control.RadioButton
        PlotBothButton               matlab.ui.control.RadioButton
        ResetDataButton              matlab.ui.control.Button
        ConstantLabel                matlab.ui.control.Label
        LoadParametersButton         matlab.ui.control.Button
        SaveParametersButton         matlab.ui.control.Button
        getesfitfit1Button           matlab.ui.control.Button
        ExportStructsButton          matlab.ui.control.Button
        DTab                         matlab.ui.container.Tab
        PausesEditFieldLabel         matlab.ui.control.Label
        PausesEditField              matlab.ui.control.NumericEditField
        plambdanPreAvgdeltadirEditFieldLabel  matlab.ui.control.Label
        ewrlsoptions                 matlab.ui.control.EditField
        ewrlsButton                  matlab.ui.control.Button
        ewrlsBusyLabel               matlab.ui.control.Label
        ApplyewrlsButton             matlab.ui.control.Button
        RevertewrlsButton            matlab.ui.control.Button
        PlayButton                   matlab.ui.control.StateButton
        CompareewrlstooriginalLabel  matlab.ui.control.Label
        ewrlsSpinner                 matlab.ui.control.Spinner
        ewrlsedLabel                 matlab.ui.control.Label
        totalslicesLabel             matlab.ui.control.Label
        mainLabel                    matlab.ui.control.Label
        VaryTabGroup                 matlab.ui.container.TabGroup
        Vary1Tab                     matlab.ui.container.Tab
        Vary1Switch                  matlab.ui.control.Switch
        Vary1Edit                    matlab.ui.control.EditField
        Vary1ListBox                 matlab.ui.control.ListBox
        Vary1RemoveButton            matlab.ui.control.Button
        Vary2Tab                     matlab.ui.container.Tab
        Vary2Switch                  matlab.ui.control.Switch
        Vary2Edit                    matlab.ui.control.EditField
        Vary2ListBox                 matlab.ui.control.ListBox
        Vary2RemoveButton            matlab.ui.control.Button
        Vary3Tab                     matlab.ui.container.Tab
        Vary3Switch                  matlab.ui.control.Switch
        Vary3Edit                    matlab.ui.control.EditField
        Vary3ListBox                 matlab.ui.control.ListBox
        Vary3RemoveButton            matlab.ui.control.Button
        Vary4Tab                     matlab.ui.container.Tab
        Vary4Switch                  matlab.ui.control.Switch
        Vary4Edit                    matlab.ui.control.EditField
        Vary4ListBox                 matlab.ui.control.ListBox
        Vary4RemoveButton            matlab.ui.control.Button
        FitOptTab                    matlab.ui.container.Tab
        FitOptEdit                   matlab.ui.control.EditField
        MethodDropDownLabel          matlab.ui.control.Label
        FitMethodDropDown            matlab.ui.control.DropDown
        TargetDropDownLabel          matlab.ui.control.Label
        TargetDropDown               matlab.ui.control.DropDown
        ScalingDropDownLabel         matlab.ui.control.Label
        FitScalingDropDown           matlab.ui.control.DropDown
        DDataButtonGroup             matlab.ui.container.ButtonGroup
        AllParallelButton            matlab.ui.control.RadioButton
        AllSequentialButton          matlab.ui.control.RadioButton
        OnlyCurrentButton            matlab.ui.control.RadioButton
        FitOptListBox                matlab.ui.control.ListBox
        SimulateButton               matlab.ui.control.Button
        Spinner                      matlab.ui.control.Spinner
        SliceLabel                   matlab.ui.control.Label
        DataFilenameLabel            matlab.ui.control.Label
        gvaluesButton                matlab.ui.control.StateButton
        SimTimeLabel                 matlab.ui.control.Label
        FittingSliceLabel            matlab.ui.control.Label
        StopFittingButton            matlab.ui.control.StateButton
        esfitButton                  matlab.ui.control.Button
        FunctionDropDown             matlab.ui.control.DropDown
        SimMethodDropDown            matlab.ui.control.DropDown
        FitGUICheckBox               matlab.ui.control.CheckBox
        RescaleDropDown              matlab.ui.control.DropDown
        RescaleLabel                 matlab.ui.control.Label
        v33Label                     matlab.ui.control.Label
        RealCheckBox                 matlab.ui.control.CheckBox
        ImaginaryCheckBox            matlab.ui.control.CheckBox
    end

    
    properties (Access = public)
        currslice=1;
        data=[];
        dataexists='n';
        simexists='n';
        whchS=[];
        fitexists='n';
        params=[];
        
        tags={};
        xlab='Field (mT)';
        inp={};
        
        SysStr={'function=','method=','weight=','S=','g=','gFrame=','lw=','lwpp=',...
            'Nucs=','n=','A_=','A=','AFrame=','tcorr=','logtcorr=','Diff=','logDiff=','DiffFrame=',...
            'lambda=','Potential=','Exchange=','J=','ee=','eeFrame=','D=','DFrame=','DStrain=',...
            'DStrainCorr=','HStrain=','gStrain=','Abund=','AStrain=','gAStrainCorr=',...
            'dvec=','eeD=','ee2=','Q=','QFrame=','nn=','nnFrame=',...
            'aF=','B0=','B2=','B4=','B6=','B8=','B10=','B12=','L=','CF0=','CF2=','CF4=','CF6=',...
            'CF8=','CF10=','CF12=','orf=','soc=','Ham='};
        ExpStr={'mwFreq=','Range=','nPoints=','ModAmp=','Mode=','CenterSweep=','Temperature=',...
            'mwCenterSweep=','mwRange=','Harmonic=','mwPhase=','mwPolarization=',...
            'CrystalOrientation=','CrystalSymmetry=','MolFrame=','Ordering='};
        StaticExpStr={'Q Value=','Power (mW)=','Power (dB)=','Reciever Gain=',...
            'Time Constant=','Conversion Time=','Averages=','Temperature='};
        OptStr={'Output=','Verbosity=','LLKM=','LLMK=','evenK=','highField=','pImax=','GridSize=','nKnots=','Transitions=',...
            'Sites=','Threshold=','Symmetry=','Enhancement=','Intensity=','Freq2Field=',...
            'IsoCutoff=','HybridCoreNuclei=','PostConvNucs=','AccumMethod='};
        rangeLH=[];
        FitOptStr={'FitRange=','maxTime=','PrintLevel=','OutArg=','RandomStart=','TolEdgeLength=',...
            'TolFun=','SimplexPars=','delta=','TolStep=','lambda=','nTrials=','PopulationSize=',...
            'maxGenerations=','GridSize=','nParticles=','SwarmParams=',};
        whchV=[];
        
        dlw=1.1;
        slw=1.5;
        plw=1.5;
        
        dclr=[0.00,0.45,0.74];
        sclr={[0.64,0.08,0.18],[0.47,0.67,0.19],[0.49,0.18,0.56],[0.93,0.69,0.13],'r'};
        pclr={[0.85,0.33,0.10],[0.49,0.18,0.56]};
        
    end
    
    
    methods (Access = private)
        
        function [B0,spec,pars,scnD] = opendata(~,name,exten)
            switch exten
                case '.DTA'
                    [B,spc,pars]=eprload(name);
                    typ='Bruker DTA';
                case '.DSC'
                    [B,spc,pars]=eprload(name);
                    typ='Bruker DTA';
                case '.YGF'
                    name=replace(name,'YGF','DTA');
                    [B,spc,pars]=eprload(name);
                    typ='Bruker DTA';
                case '.spc'
                    [B,spcH,pars]=eprload(name);
                    if iscell(B)
                        B{1,2}=B{1,2}';
                        spc=spcH;
                    else
                        B=B';
                        spc=spcH';
                    end
                    typ='Bruker spc';
                case '.par'
                    [B,spcH,pars]=eprload(name);
                    if iscell(B)
                        B{1,2}=B{1,2}';
                        spc=spcH;
                    else
                        B=B';
                        spc=spcH';
                    end
                    typ='Bruker spc';
                case '.ESR'
                    [B,spc,pars]=eprload(name);
                    typ='ActiveSpectrum';
                case '.xml'
                    [BNaN,spcNaN,pars]=eprload(name);
                    BNaN=BNaN.*10;
                    m=1;
                    for i=1:length(BNaN)
                        if ~isnan(BNaN(i)) && BNaN(i)~=0
                            Bm(m,1)=BNaN(i);
                            spcm(m,1)=spcNaN(i);
                            m=m+1;
                        else
                        end
                    end
                    B=Bm;
                    spc=spcm;
                    typ='magnettech';
                otherwise
                    original=readmatrix(name);
                    B=original(:,1);
                    spc=original(:,2:end);
                    typ='ascii';
            end
            
            if length(spc(1,:))>1
                switch typ
                    case 'ascii'
                        B_(:,1)=B;
                        scnD1=1:1:length(spc(1,:));
                    otherwise
                        B_(:,1)=B{1,1};
                        scnD1=B{1,2};
                end
                k=1;
                %assignin("base",'spc',spc)
                for i=1:length(spc(1,:))
                    if sum(abs(spc(:,i)))~=0
                        spk(:,k)=spc(:,i);
                        scnD(k)=scnD1(i);
                        k=k+1;
                    else
                    end
                end
                %assignin("base",'spk',spk)
            else
                B_=B;
                spk=spc;
                scnD=[];
            end
            
            B1=B_;
            for i=2:length(B_)+1
                if i>=length(B_)
                    break
                else
                    if B_(i)==B_(i-1) || B_(i)==B_(i+1)
                        B1(i,1)=(B_(i-1)+B_(i+1))/2;
                    else
                        B1(i,1)=B_(i);
                    end
                end
            end
            
            if max(B1)<1000
                B0=B1;
            else
                B0=B1/10;
            end
            
            if length(spk(1,:))>1
                for i=1:length(spk(1,:))
                    [spec(:,i),~]=basecorr(spk(:,i),1,0);
                end
            else
                [spec,~]=basecorr(spk,1,0);
            end
            
            try
                pars.MFQ=num2str(pars.MWFQ/1e9);
            catch
                try
                    pars.MFQ=num2str(pars.MF);
                catch
                    try
                        pars.MFQ=num2str(pars.MwFreq);
                    catch
                        pars.MFQ='';
                    end
                end
            end
            
            try
                [modA_,~]=strtok(pars.ModAmp,' G');
                modA=str2double(modA_)/10;
                pars.modA=num2str(modA);
            catch
                try
                    [modA,~]=strtok(pars.Modulation,' G');
                    pars.modA=num2str(modA/10);
                catch
                    try
                        pars.modA=num2str(pars.RMA/10);
                    catch
                        pars.modA='';
                    end
                end
            end
            
            try
                [ct,~]=strtok(pars.ConvTime,' ms');
                pars.CT=num2str(ct);
            catch
                try
                    pars.CT=num2str(pars.RCT);
                catch
                    pars.CT='';
                end
            end
            
            try
                [tc,~]=strtok(pars.TimeConst,' ms');
                pars.TC=num2str(tc);
            catch
                try
                    pars.TC=num2str(pars.RTC);
                catch
                    pars.TC='';
                end
            end
            
            try
                [rg,~]=strtok(pars.Gain,' dB');
                pars.RG=num2str(rg);
            catch
                try
                    pars.RG=num2str(pars.RRG);
                catch
                    pars.RG='';
                end
            end
            
            try
                [pow,~]=strtok(pars.Power,' mW');
                pars.POW=num2str(pow);
            catch
                try
                    pars.POW=num2str(pars.MP);
                catch
                    pars.POW='';
                end
            end
            
            try
                [attn,~]=strtok(pars.PowerAtten,' dB');
                pars.ATTN=num2str(attn);
            catch
                try
                    pars.ATTN=num2str(pars.MPD);
                catch
                    pars.ATTN='';
                end
            end
            
            try
                pars.QVal=num2str(pars.QValue);
            catch
                pars.QVal='';
            end
            
            try
                pars.Harm=num2str(pars.Harmonic);
            catch
                pars.Harm='';
            end
            
            try
                pars.nScans=num2str(pars.AVGS);
            catch
                pars.nScans='';
            end
            
            try
                [temper,~]=strtok(pars.Temperature,' K');
                pars.Temp = num2str(temper);
            catch
                pars.Temp = '';
            end
            
        end
        
        function [sys,fncmeth_f,fncmeth_m]=getsyss(~,syslst)
            
            l=1;
            SysValues=[];
            for i=1:numel(syslst)
                syslst_f=strsplit(syslst{i},'=');
                syslst_f{2} = erase(syslst_f{2},"'");
                switch syslst_f{1}
                    case 'function'
                        fncmeth_f=syslst_f{2};
                    case 'method'
                        fncmeth_m=syslst_f{2};
                    otherwise
                        if ~isempty(syslst_f{2})
                            SysString{l}=syslst_f{1};
                            SysValues{l}=str2num(syslst_f{2});
                            if isempty(SysValues{l})
                                SysValues{l}=syslst_f{2};
                            end
                            l=l+1;
                        else
                        end
                end
            end
            if isempty(SysValues)
                sys=[];
            else
                sys=cell2struct(SysValues,SysString,2);
            end
            
        end
        
        function [Exp] = getexps(~,explst)
            
            h=1;
            ExpValues = [];
            for i=1:numel(explst)
                explst_f=strsplit(explst{i},'=');
                explst_f{2}=erase(explst_f{2},"'");
                if ~isempty(explst_f{2})
                    ExpString{h}=explst_f{1};
                    ExpValues{h}=str2num(explst_f{2});
                    if isempty(ExpValues{h})
                        ExpValues{h}=explst_f{2};
                    end
                    h=h+1;
                end
            end
            if isempty(ExpValues)
                Exp = [];
            else
                Exp=cell2struct(ExpValues,ExpString,2);
            end
            
        end
        
        function [Opt] = getopts(~,optlst)
            
            k=1;
            OptValues=[];
            for i=1:numel(optlst)
                optlst_f=strsplit(optlst{i},'=');
                optlst_f{2}=erase(optlst_f{2},"'");
                if ~isempty(optlst_f{2})
                    OptString{k}=optlst_f{1};
                    OptValues{k}=str2num(optlst_f{2});
                    if isempty(OptValues{k})
                        OptValues{k}=optlst_f{2};
                    end
                    k=k+1;
                end
            end
            if isempty(OptValues)
                Opt=[];
            else
                Opt=cell2struct(OptValues,OptString,2);
            end
            
        end
        
        function [syslst_new]=sysadd(~,sysstr,syslst)
            
            sysparam=strsplit(sysstr,'=');
            syslst_new=syslst;
            
            for i=1:numel(syslst)
                syslstf=strsplit(syslst{i},'=');
                switch syslstf{1}
                    case sysparam{1}
                        syslst_new{i}=[sysparam{1},'=',sysparam{2}];
                    otherwise
                        syslst_new{i}=syslst{i};
                end
            end
            
        end
        
        function syslst_new = getsysstructs(~,syslst,sys_strct,sys_params)
            
            syslst_new=syslst;
            for i=1:numel(sys_params)
                for j=1:numel(syslst)
                    syslstf=strsplit(syslst{j},'=');
                    switch syslstf{1}
                        case sys_params{i}
                            if contains(sys_params{i},'function')==1 || contains(sys_params{i},'method')==1 || contains(sys_params{i},'Nucs')==1
                                syslst_new{j}=[sys_params{i},'=',sys_strct.(sys_params{i})];
                            else
                                sys_vals{1}=sys_params{i};
                                sys_val=sys_strct.(sys_params{i});
                                if length(sys_val(1,:))==1 && length(sys_val(:,1))==1
                                    sys_vals{2}=num2str(sys_val(1));
                                    
                                elseif length(sys_val(1,:))==2 && length(sys_val(:,1))==1
                                    sys_vals{2}=['[',num2str(sys_val(1)),' ',num2str(sys_val(2)),']'];
                                    
                                elseif length(sys_val(1,:))==1 && length(sys_val(:,1))==2
                                    sys_vals{2}=['[',num2str(sys_val(1)),';',num2str(sys_val(2)),']'];
                                    
                                elseif length(sys_val(1,:))==2 && length(sys_val(:,1))==2
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),';',num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),']'];
                                    
                                elseif length(sys_val(1,:))==3 && length(sys_val(:,1))==1
                                    sys_vals{2}=['[',num2str(sys_val(1)),' ',num2str(sys_val(2)),' ',num2str(sys_val(3)),']'];
                                    
                                elseif length(sys_val(1,:))==1 && length(sys_val(:,1))==3
                                    sys_vals{2}=['[',num2str(sys_val(1)),';',num2str(sys_val(2)),';',num2str(sys_val(3)),']'];
                                    
                                elseif length(sys_val(1,:))==3 && length(sys_val(:,1))==2
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),';',...
                                        num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',num2str(sys_val(2,3)),']'];
                                    
                                elseif length(sys_val(1,:))==2 && length(sys_val(:,1))==3
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),'; ',num2str(sys_val(2,1)),' ',...
                                        num2str(sys_val(2,2)),'; ',num2str(sys_val(3,1)),' ',num2str(sys_val(3,2)),']'];
                                    
                                elseif length(sys_val(1,:))==3 && length(sys_val(:,1))==3
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),';',...
                                        num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',num2str(sys_val(2,3)),';',...
                                        num2str(sys_val(3,1)),' ',num2str(sys_val(3,2)),' ',num2str(sys_val(3,3)),']'];
                                    
                                elseif length(sys_val(1,:))==5 && length(sys_val(:,1))==2
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),...
                                        ' ',num2str(sys_val(1,3)),' ',num2str(sys_val(1,4)),' ',num2str(sys_val(1,5)),...
                                        ';',num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',...
                                        num2str(sys_val(2,3)),' ',num2str(sys_val(2,4)),' ',num2str(sys_val(2,5)),']'];
                                    
                                elseif length(sys_val(1,:))==5 && length(sys_val(:,1))==1
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),' ',...
                                        num2str(sys_val(1,4)),' ',num2str(sys_val(1,5)),']'];
                                    
                                elseif length(sys_val(1,:))==6 && length(sys_val(:,1))==2
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),...
                                        ' ',num2str(sys_val(1,3)),' ',num2str(sys_val(1,4)),' ',num2str(sys_val(1,5)),...
                                        ' ',num2str(sys_val(1,6)),';',num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',...
                                        num2str(sys_val(2,3)),' ',num2str(sys_val(2,4)),' ',num2str(sys_val(2,5)),' ',num2str(sys_val(2,6)),']'];
                                    
                                elseif length(sys_val(1,:))==6 && length(sys_val(:,1))==1
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),' ',...
                                        num2str(sys_val(1,4)),' ',num2str(sys_val(1,5)),' ',num2str(sys_val(1,6)),']'];
                                    
                                elseif length(sys_val(1,:))==9 && length(sys_val(:,1))==1
                                    sys_vals{2}=['[',num2str(sys_val(1)),' ',num2str(sys_val(2)),' ',num2str(sys_val(3)),' ',...
                                        num2str(sys_val(4)),' ',num2str(sys_val(5)),' ',num2str(sys_val(6)),' ',...
                                        num2str(sys_val(7)),' ',num2str(sys_val(8)),' ',num2str(sys_val(9)),']'];
                                    
                                elseif length(sys_val(1,:))==3 && length(sys_val(:,1))==6
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),';',...
                                        num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',num2str(sys_val(2,3)),';',...
                                        num2str(sys_val(3,1)),' ',num2str(sys_val(3,2)),' ',num2str(sys_val(3,3)),']; [',...
                                        num2str(sys_val(4,1)),' ',num2str(sys_val(4,2)),' ',num2str(sys_val(4,3)),';',...
                                        num2str(sys_val(5,1)),' ',num2str(sys_val(5,2)),' ',num2str(sys_val(5,3)),';',...
                                        num2str(sys_val(6,1)),' ',num2str(sys_val(6,2)),' ',num2str(sys_val(6,3)),']'];
                                    
                                elseif length(sys_val(1,:))==6 && length(sys_val(:,1))==3
                                    sys_vals{2}=['[',num2str(sys_val(1,1)),' ',num2str(sys_val(1,2)),' ',num2str(sys_val(1,3)),';',...
                                        num2str(sys_val(2,1)),' ',num2str(sys_val(2,2)),' ',num2str(sys_val(2,3)),';',...
                                        num2str(sys_val(3,1)),' ',num2str(sys_val(3,2)),' ',num2str(sys_val(3,3)),'], [',...
                                        num2str(sys_val(1,4)),' ',num2str(sys_val(1,5)),' ',num2str(sys_val(1,6)),';',...
                                        num2str(sys_val(2,4)),' ',num2str(sys_val(2,5)),' ',num2str(sys_val(2,6)),';',...
                                        num2str(sys_val(3,4)),' ',num2str(sys_val(3,5)),' ',num2str(sys_val(3,6)),']'];
                                end
                                syslst_fil=[sys_vals(1),'=',sys_vals(2)];
                                syslst_new{j}=[syslst_fil{1},syslst_fil{2},syslst_fil{3}];
                            end
                        otherwise
                            syslst_new{j}=syslst{j};
                    end
                end
                syslst=syslst_new;
            end
            
        end
        
        function [vary]=getvarys(~,varypars,t)
            
            j=1;
            for i=1:t
                varypar=strsplit(varypars{i},'=');
                VaryString{j}=varypar{1};
                VaryValues{j}=str2num(varypar{2});
                j=j+1;
            end
            
            vary=cell2struct(VaryValues,VaryString,2);
            
        end
        function [FitOpt, fitrange] = getfitopts(~,fitoptlst)
            
            l=1;
            fitOptValues=[];
            for i=1:numel(fitoptlst)
                fitoptlst_f=strsplit(fitoptlst{i},'=');
                fitoptlst_f{2}=erase(fitoptlst_f{2},"'");
                switch fitoptlst_f{1}
                    case 'FitRange'
                        fitrange=str2num(fitoptlst_f{2});
                    otherwise
                        if ~isempty(fitoptlst_f{2})
                            fitOptString{l}=fitoptlst_f{1};
                            fitOptValues{l}=str2num(fitoptlst_f{2});
                            if isempty(fitOptValues{l})
                                fitOptValues{l}=fitoptlst_f{2};
                            end
                            l=l+1;
                        end
                end
            end
            if isempty(fitOptValues)
                FitOpt=[];
            else
                FitOpt=cell2struct(fitOptValues,fitOptString,2);
            end
            
        end
        
        function [varylst_new]=varyadd(~,varystr,varylst1)
            
            if isempty(varylst1)
                varylst_new{1}=varystr;
            else
                varylst_new=varylst1;
                varylst=varylst1;
                
                varyparam=strsplit(varystr,'=');
                if ~isempty(varylst1)
                    vn=numel(varylst);
                    n=0;
                    for i=1:vn
                        varylstf=strsplit(varylst{i},'=');
                        switch varylstf{1}
                            case varyparam{1}
                                varylst_new{i}=[varylstf{1},'=',varyparam{2}];
                            otherwise
                                varylst_new{i}=varylst{i};
                                n=n+1;
                        end
                    end
                    if n==vn
                        varylst_new{vn+1}=[varyparam{1},'=',varyparam{2}];
                    end
                end
            end
            
        end
        
        function [parlst_new]=paradd(~,parstr,parlst)
            
            parparam=strsplit(parstr,'=');
            parlst_new=parlst;
            
            for i=1:numel(parlst)
                parlstf=strsplit(parlst{i},'=');
                switch parlstf{1}
                    case parparam{1}
                        parlst_new{i}=[parparam{1},'=',parparam{2}];
                    otherwise
                        parlst_new{i}=parlst{i};
                end
            end
            
        end
        
        function [simX,sim]=custom_twomethods(~,Sys,Exp,options)
            
            l=length(Sys);
            sim=zeros(Exp.nPoints,l);
            if l==1
                Opt = options;
                switch options.fnc{1}
                    case 'garlic'
                        Opt.Method=options.meth{1};
                        [simX(:,1),sim(:,1)] = garlic(Sys,Exp,Opt);
                    case 'chili'
                        Opt.LiouvMethod=options.meth{1};
                        [simX(:,1),sim(:,1)] = chili(Sys,Exp,Opt);
                    case 'pepper'
                        Opt.Method=options.meth{1};
                        [simX(:,1),sim(:,1)] = pepper(Sys,Exp,Opt);
                end
                %simout(:,1)=sim(Opt.rangeL:Opt.rangeH);
            else
                for i=1:l
                    Opt = options;
                    switch options.fnc{i}
                        case 'garlic'
                            Opt.Method=options.meth{i};
                            [simX(:,1),sim(:,i)] = garlic(Sys{i},Exp,Opt);
                        case 'chili'
                            Opt.LiouvMethod=options.meth{i};
                            [simX(:,1),sim(:,i)] = chili(Sys{i},Exp,Opt);
                        case 'pepper'
                            Opt.Method=options.meth{i};
                            [simX(:,1),sim(:,i)] = pepper(Sys{i},Exp,Opt);
                    end
                    %simout(:,i)=sim(Opt.rangeL:Opt.rangeH);
                end
            end
            
        end
        
        function simout=custom_fittwomethods(~,Sys,Exp,options)
            
            l=length(Sys);
            sim=zeros(Exp.nPoints,l);
            if l==1
                Opt = options;
                switch options.fnc{1}
                    case 'garlic'
                        Opt.Method=options.meth{1};
                        [~,sim] = garlic(Sys,Exp,Opt);
                    case 'chili'
                        Opt.LiouvMethod=options.meth{1};
                        [~,sim] = chili(Sys,Exp,Opt);
                    case 'pepper'
                        Opt.Method=options.meth{1};
                        [~,sim] = pepper(Sys,Exp,Opt);
                end
                simspc=sim;
            else
                for i=1:l
                    Opt = options;
                    switch options.fnc{i}
                        case 'garlic'
                            Opt.Method=options.meth{i};
                            [~,sim(:,i)] = garlic(Sys{i},Exp,Opt);
                        case 'chili'
                            Opt.LiouvMethod=options.meth{i};
                            [~,sim(:,i)] = chili(Sys{i},Exp,Opt);
                        case 'pepper'
                            Opt.Method=options.meth{i};
                            [~,sim(:,i)] = pepper(Sys{i},Exp,Opt);
                    end
                    %sim(:,i)=system.weight.*(sim(:,i)./max(sim(:,i)));
                end
                simspc=sum(sim,2);
            end
            simout=zeros(length(simspc),1);
            simout(options.rangeLH(1):options.rangeLH(2))=...
                simspc(options.rangeLH(1):options.rangeLH(2));
        end
        
        function simout=custom_fit2D(~,Sys,Exp,options)
            
            l=length(Sys);
            simout=zeros(1,Exp.nP*Exp.nSpcs);
            Exp.nPoints=Exp.nP;
            for i=1:Exp.nSpcs
                k=((i-1)*Exp.nP)+1;
                kk=((i)*Exp.nP);
                sim=zeros(Exp.nP,l);
                if l==1
                    Opt = options;
                    switch options.fnc{1}
                        case 'garlic'
                            Opt.Method=options.meth{1};
                            [~,sim] = garlic(Sys,Exp,Opt);
                        case 'chili'
                            Opt.LiouvMethod=options.meth{1};
                            [~,sim] = chili(Sys,Exp,Opt);
                        case 'pepper'
                            Opt.Method=options.meth{1};
                            [~,sim] = pepper(Sys,Exp,Opt);
                    end
                    simSpc=sim./max(sim);
                else
                    for j=1:l
                        Opt = options;
                        switch options.fnc{j}
                            case 'garlic'
                                Opt.Method=options.meth{j};
                                [~,sim(:,j)] = garlic(Sys{j},Exp,Opt);
                            case 'chili'
                                Opt.LiouvMethod=options.meth{j};
                                [~,sim(:,j)] = chili(Sys{j},Exp,Opt);
                            case 'pepper'
                                Opt.Method=options.meth{i};
                                [~,sim(:,j)] = pepper(Sys{j},Exp,Opt);
                        end
                        %     sim(:,i)=str2num(Opt.wght{i}).*(sim(:,i)./max(sim(:,i)));
                    end
                    simspc=sum(sim,2);
                    simSpc=simspc./max(simspc);
                end
                sim_spc=zeros(length(simSpc),1);
                sim_spc(options.rangeLH(1):options.rangeLH(2))=...
                    simSpc(options.rangeLH(1):options.rangeLH(2));
                simout(1,k:kk)=sim_spc;
            end
        end
        
        %{
        function simout=custom_fitALL(~,Sys,Exp,Opt)
            
            Exp.nPoints = Exp.nP;
            for i=1:Exp.nSpcs
                k=((i-1)*Exp.nP)+1;
                kk=((i)*Exp.nP);
                switch Exp.func
                    case 'garlic'
                        [~,simSpc]=garlic(Sys,Exp,Opt);
                    case 'chili'
                        [~,simSpc]=chili(Sys,Exp,Opt);
                    case 'pepper'
                        [~,simSpc]=pepper(Sys,Exp,Opt);
                    case 'salt'
                        [~,simSpc]=salt(Sys,Exp,Opt);
                    case 'saffron'
                        [~,simSpc]=saffron(Sys,Exp,Opt);
                end
                simout(1,k:kk)=simSpc;
            end
            
        end
        %}
        
        function showplots(app,selbutton,whch,onORoff)
            
            switch app.simexists
                case 'y'
                    l=length(whch);
                    if l~=0
                        k=1;
                        for i=1:4
                            switch onORoff{i}
                                case 'On'
                                    inti(k)=max(app.data.simspc(:,i));
                                    simspci(:,k)=app.data.simspc(:,i);
                                    k=k+1;
                            end
                        end
                        totsimspc_=sum(simspci,2);
                        intsim=max(totsimspc_);
                        totsimspc=totsimspc_./intsim;
                        for i=1:length(simspci(1,:))
                            simspci(:,i)=(simspci(:,i)./max(app.data.simspc(:,i))).*(inti(i)/intsim);
                        end
                    end
                case 'n'
                    selbutton='Data Only';
                    set(app.SimOnlyButton,'Value',0)
                    set(app.DataOnlyButton,'Value',1)
                    set(app.PlotBothButton,'Value',0)
            end
            switch app.dataexists
                case 'y'
                    dataspc=app.data.spc(:,app.currslice);
                    switch selbutton
                        case 'Data Only'
                            plot(app.UIAxes,app.data.B,dataspc,'linewidth',app.dlw,'color',app.dclr);
                            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
                            ylim(app.UIAxes,[min(dataspc) max(dataspc)].*1.05);
                            legend(app.UIAxes,'Data','location','northeast')
                        case 'Sim Only'
                            if l~=0
                                for i=1:l
                                    plot(app.UIAxes,app.data.simB, simspci(:,i),'linewidth',app.slw,'color',app.sclr{i})
                                    lg{i}=['SYS',num2str(whch(i))];
                                    hold(app.UIAxes,'on');
                                end
                                if l>1
                                    plot(app.UIAxes,app.data.simB,totsimspc,'linewidth',app.slw,'color',app.sclr{5});
                                    lg{l+1}='Sim';
                                end
                                xlim(app.UIAxes,[min(app.data.simB) max(app.data.simB)]);
                                ylim(app.UIAxes,[min(totsimspc) max(totsimspc)].*1.05);
                                legend(app.UIAxes,lg,'location','northeast')
                                hold(app.UIAxes,'off');
                            end
                        case 'Both'
                            lg{1}='Data';
                            plot(app.UIAxes,app.data.B,dataspc,'linewidth',app.dlw,'color',app.dclr);
                            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
                            ylim(app.UIAxes,[min(dataspc) max(dataspc)].*1.05);
                            legend(app.UIAxes,lg,'location','northeast')
                            if l~=0
                                switch app.RescaleDropDown.Value
                                    case 'none'
                                        totsimspc=totsimspc_;
                                    otherwise
                                        try
                                            totsimspc=rescaledata(totsimspc_,dataspc, app.RescaleDropDown.Value);
                                        catch
                                            totsimspc=rescale(totsimspc_,dataspc, app.RescaleDropDown.Value);
                                        end
                                end
                                hold(app.UIAxes,'on');
                                %{
                                  for i=1:l
                                     plot(app.UIAxes,app.data.simB, simspci(:,i),'linewidth',app.slw,'color',app.sclr{i})
                                     lg{1+i}=['SYS',num2str(whch(i))];
                                  end
                                  if l>1
                                %}
                                plot(app.UIAxes,app.data.B,totsimspc,'linewidth',app.slw,'color',app.sclr{5});
                                lg{2}='Sim';
                                %lg{l+2}='Sim';
                                %                            end
                                legend(app.UIAxes,lg,'location','northeast')
                                hold(app.UIAxes,'off');
                                sim_rmsd = sqrt(mean((dataspc-totsimspc).^2));
                                text(app.UIAxes,0.85,0.05,['RMSD = ',num2str(sim_rmsd)],'Units','normalized')
                            end
                    end
                case 'n'
                    set(app.SimOnlyButton,'Value',1)
                    set(app.DataOnlyButton,'Value',0)
                    set(app.PlotBothButton,'Value',0)
                    if l~=0
                        for i=1:l
                            plot(app.UIAxes,app.data.simB, simspci(:,i),'linewidth',app.slw,'color',app.sclr{i})
                            lg{i}=['SYS',num2str(whch(i))];
                            hold(app.UIAxes,'on');
                        end
                        if l>1
                            plot(app.UIAxes,app.data.simB,totsimspc,'linewidth',app.slw,'color',app.sclr{5});
                            lg{l+1}='Sim';
                        end
                        xlim(app.UIAxes,[min(app.data.simB) max(app.data.simB)]);
                        ylim(app.UIAxes,[min(totsimspc) max(totsimspc)].*1.05);
                        legend(app.UIAxes,lg,'location','northeast')
                        hold(app.UIAxes,'off');
                    end
            end
            xlabel(app.UIAxes,app.xlab);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
        end
        
        function togglesysplts(app,onORoff,selbutton)
            
            switch app.simexists
                case 'y'
                    whchSplt=app.whchS;
                    if contains(num2str(app.whchS),'1')
                        onORoff{1}=get(app.Sys1Switch,'Value');
                        switch onORoff{1}
                            case 'Off'
                                whchSplt=str2num(erase(num2str(whchSplt),'1'));
                        end
                    end
                    if contains(num2str(app.whchS),'2')
                        onORoff{2}=get(app.Sys2Switch,'Value');
                        switch onORoff{2}
                            case 'Off'
                                whchSplt=str2num(erase(num2str(whchSplt),'2'));
                        end
                    end
                    if contains(num2str(app.whchS),'3')
                        onORoff{3}=get(app.Sys3Switch,'Value');
                        switch onORoff{3}
                            case 'Off'
                                whchSplt=str2num(erase(num2str(whchSplt),'3'));
                        end
                    end
                    if contains(num2str(app.whchS),'4')
                        onORoff{4}=get(app.Sys4Switch,'Value');
                        switch onORoff{4}
                            case 'Off'
                                whchSplt=str2num(erase(num2str(whchSplt),'4'));
                        end
                    end
                    if isempty(whchSplt)
                        switch app.dataexists
                            case 'y'
                                selbutton='Data Only';
                            case 'n'
                                return
                        end
                    else
                        switch app.dataexists
                            case 'n'
                                selbutton='Sim Only';
                        end
                    end
                case 'n'
                    whchSplt=[];
                    switch app.dataexists
                        case 'y'
                            selbutton='Data Only';
                        case 'n'
                            return
                    end
            end
            showplots(app,selbutton,whchSplt,onORoff)
            
        end
        
        function sliceplts(app)
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            switch app.tags.lastprocess
                case 'baseline'
                    plot(app.UIAxes,app.data.B,app.data.subspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
                    plot(app.UIAxes,app.data.B,app.data.fitln(:,app.currslice),'linewidth',app.plw,'color',app.pclr{2});
                    ylim(app.UIAxes,[min(app.data.subspc(:,app.currslice)) max(app.data.subspc(:,app.currslice))].*1.05);
                    legend(app.UIAxes,'Original','Baseline','Corrected','location','northeast');
                    
                case 'smoothing'
                    plot(app.UIAxes,app.data.B,app.data.fspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
                    legend(app.UIAxes,'Original','Filtered','location','northeast');
                    
                case 'integral'
                    set(app.DIEditField,'Value',num2str(app.data.dblintgrl(app.currslice)))
                    plt_fac_dbl = (max(app.data.spc(:,app.currslice))/max(app.data.dblintgrted(:,app.currslice)));
                    plt_fac_sngl = (max(app.data.spc(:,app.currslice))/max(app.data.intgrted(:,app.currslice)));
                    switch app.tags.ShowInt
                        case 'Single'
                            plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                            legend(app.UIAxes,'Data','First Integral','location','northwest');
                        case 'Double'
                            plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{1});
                            legend(app.UIAxes,'Data','Double Integral','location','northwest');
                        case 'Both'
                            plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                            plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{2});
                            legend(app.UIAxes,'Data','First Integral','Double Integral','location','northwest');
                    end
                    
                case 'interpolated'
                    plot(app.UIAxes,app.data.B,app.data.interpspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
                    legend(app.UIAxes,'Original','Interpolated','location','northeast');
                    
                case 'noised'
                    plot(app.UIAxes,app.data.B,app.data.noisespc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
                    legend(app.UIAxes,'Original','Noised','location','northeast');
                    
                case 'simulation'
                    totsimspc_=sum(app.data.simspc,2);
                    switch app.RescaleDropDown.Value
                        case 'none'
                            totsimspc=totsimspc_;
                        otherwise
                            try
                                totsimspc=rescaledata(totsimspc_,app.data.spc(:,app.currslice),app.RescaleDropDown.Value);
                            catch
                                totsimspc=rescale(totsimspc_,app.data.spc(:,app.currslice),app.RescaleDropDown.Value);
                            end
                    end
                    plot(app.UIAxes,app.data.B,totsimspc,'linewidth',app.slw,'color',app.sclr{5});
                    legend(app.UIAxes,'Data','Simulation','location','northeast');
                    sim_rmsd = sqrt(mean((app.data.spc(:,app.currslice)-totsimspc).^2));
                    text(app.UIAxes,0.85,0.05,['RMSD = ',num2str(sim_rmsd)],'Units','normalized')
            end
            hold(app.UIAxes,'off');
            
            switch app.tags.secondmomented
                case 'y'
                    set(app.SecondMomentEdit,'Value',num2str(app.data.second_moment(:,app.currslice)))
            end
            
            switch app.tags.spincounted
                case 'y'
                    set(app.SpinCEditField,'Value',num2str(app.data.spinC(:,app.currslice)))
            end
            
        end
        
        function slide_freq(app,slid_val)
            
            switch app.tags.gvaled
                case 'y'
                    app.data.B=app.data.ogB;
                    set(app.UIAxes,'xdir','normal');
                    set(app.gvaluesButton,'Text','g values','Value',0)
                    app.tags.gvaled='n';
            end
            
            gaxis = (((6.6261e-34)*app.params.simFq*1e9)./((9.2740e-24).*app.data.B)).*1e3;
            
            B = (((6.6261e-34)*slid_val*1e9)./((9.2740e-24).*gaxis)).*1e3;
            
            totsimspc=sum(app.data.simspc,2);
            switch app.RescaleDropDown.Value
                case 'none'
                    simspc=totsimspc./max(totsimspc);
                otherwise
                    try
                        simspc=rescaledata(totsimspc,app.data.spc(:,app.currslice), app.RescaleDropDown.Value);
                    catch
                        simspc=rescale(totsimspc,app.data.spc(:,app.currslice), app.RescaleDropDown.Value);
                    end
            end
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on')
            plot(app.UIAxes,B,simspc,'linewidth',app.slw,'color',app.sclr{5});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,'Field (mT)');
            legend(app.UIAxes,'Data','Simulation','location','northeast');
            hold(app.UIAxes,'off')
            
        end
        
        function [intgrted, dblintgrted, intgrl, dblintgrl, l1, l2] =...
                integrate_calculate(app, baseperc, spec)
            
            baseselect=get(app.BaselineCorrButtonGroup,'SelectedObject');
            spc=basecorr(spec,1,0);
            b_inc = mean(diff(app.data.B)*10);
            
            l=length(spc);
            if baseperc == 0 || baseperc == 50
                l1 = 1;
                l2 = l;
            else
                l1=round(l*baseperc/100);
                l2=l-l1;
            end
            polyorder=get(app.IntOrderSpinner,'Value');
            if isempty(polyorder)
                polyorder=0;
            end
            switch app.TrapzCheckBox.Value
                case 1
                    for i=1:length(spc(1,:))
                        frstint(:,i)=cumtrapz(b_inc,spc(:,i));
                        frstintgrl(i)=trapz(b_inc,frstint(:,i));
                    end
                case 0
                    for i=1:length(spc(1,:))
                        frstint(:,i)=cumsum(spc(:,i)).*b_inc;
                        frstintgrl(i)=sum(frstint(:,i))*b_inc;
                    end
            end
            
            for i=1:length(frstint(1,:))
                truncspc1=frstint(1:l1,i);
                truncspc2=frstint(l-l1+1:l,i);
                truncspc3=linspace(frstint(l1,i),frstint(l-l1+1,i),l-(2*l1))';
                truncspc(1:l1)=truncspc1;
                truncspc(l1+1:l-l1)=truncspc3;
                truncspc(l-l1+1:l)=truncspc2;
                switch baseselect.Text
                    case 'None'
                        intgrted(:,i)=frstint(:,i);
                    case 'Polynomial, order:'
                        [~,fitln]=basecorr(truncspc',1,polyorder);
                        intgrted(:,i)=(frstint(:,i)-fitln);
                end
                switch app.TrapzCheckBox.Value
                    case 1
                        intgrl(i)=trapz(b_inc,intgrted(:,i));
                    case 0
                        intgrl(i)=sum(intgrted(:,i))*b_inc;
                end
            end
            
            switch app.TrapzCheckBox.Value
                case 1
                    for i=1:length(intgrted(1,:))
                        dblintgrted(:,i)=cumtrapz(b_inc,intgrted(:,i));
                        dblintgrl(i)=trapz(b_inc,intgrted(l1:l2,i));
                    end
                case 0
                    for i=1:length(intgrted(1,:))
                        dblintgrted(:,i)=cumsum(intgrted(:,i))*b_inc;
                        dblintgrl(i)=sum(intgrted(l1:l2,i))*b_inc;
                    end
            end
            
        end
        
    end
    

    methods (Access = private)

        % Code that executes after component creation
        function Opening(app)
            
            scrn = get(groot,'ScreenSize');
            drawnow
            app.cwEPRapp.Position = ...
                [scrn(3)*.07 scrn(4)*.12 scrn(3)*.8333 scrn(4)*.8];
            
            
            set(app.Sys1ListBox,'Items',app.SysStr)
            set(app.Sys2ListBox,'Items',app.SysStr)
            set(app.Sys3ListBox,'Items',app.SysStr)
            set(app.Sys4ListBox,'Items',app.SysStr)
            
            set(app.ExpListBox,'Items',app.ExpStr)
            set(app.StaticExpListBox,'Items',app.StaticExpStr)
            set(app.OptListBox,'Items',app.OptStr)
            set(app.FitOptListBox,'Items',app.FitOptStr)
            
            set(app.Sys1Switch,'Value','On')
            set(app.Sys2Switch,'Value','Off')
            set(app.Sys3Switch,'Value','Off')
            set(app.Sys4Switch,'Value','Off')
            set(app.Vary1Switch,'Value','Off')
            set(app.Vary2Switch,'Value','Off')
            set(app.Vary3Switch,'Value','Off')
            set(app.Vary4Switch,'Value','Off')
            
            %set(app.DIfactorEditField,'Value','Power,Q,ModAmp,S,Temp')
            
        end

        % Button pushed function: OpenDataButton
        function Opendata(app, event)
            app.cwEPRapp.Visible = 'Off';
            try
                [name,path]=uigetfile({'*.*';'*.DTA';'*.spc';'*.ESR';'*.txt';'*.dat';'*.xml'});
                [rootnm,exten]=strtok(name,'.');
                set(app.DataFilenameLabel,'Text',rootnm)
            catch
                app.cwEPRapp.Visible = 'On';
                return
            end
            app.cwEPRapp.Visible = 'On';
            
            app.data.B=[];
            app.data.spc = [];
            app.data.orgB=[];
            app.data.orgrealspc=[];
            app.data.orgimagspc = [];
            
            [app.data.B,spc,app.params,app.data.scnD] =...
                opendata(app,strcat(path,filesep,name),exten);
            
            set(app.RealCheckBox,'Value',1)
            set(app.ImaginaryCheckBox,'Value',0)
            if ~isreal(spc)
                app.RealCheckBox.Visible = 'On';
                app.ImaginaryCheckBox.Visible = 'On';
                app.data.complexspc = spc;
                app.data.orgrealspc = real(spc);
                app.data.orgimagspc = imag(spc);
            elseif isreal(spc)
                app.RealCheckBox.Visible = 'Off';
                app.ImaginaryCheckBox.Visible = 'Off';
            end
            
            app.data.spc = real(spc);
            app.data.orgspc = app.data.spc;
            app.data.orgB=app.data.B;
            app.currslice=1;
            set(app.PointsSpinner,'Value',length(app.data.spc(:,1)))
            
            if length(app.data.spc(1,:))>1
                set(app.Spinner,'Visible','On')
                set(app.Spinner,'Limits',[1 length(app.data.spc(1,:))],'Value',1)
                set(app.SliceLabel,'Visible','On')
                set(app.totalslicesLabel,'Visible','on','Text',[num2str(length(app.data.spc(1,:))),' slices'])
            else
                set(app.Spinner,'Visible','Off','Value',1)
                set(app.SliceLabel,'Visible','Off')
                set(app.totalslicesLabel,'Visible','on','Text','1D Data')
            end
            app.tags.lastprocess='none';
            app.tags.baselined='n';
            app.tags.subtracted='n';
            app.tags.filtered='n';
            app.tags.smoothed='n';
            app.tags.integrated='n';
            set(app.DIEditField,'Value','')
            app.tags.secondmomented='n';
            set(app.SecondMomentEdit,'Value','')
            app.tags.spincounted = 'n';
            set(app.SpinCEditField,'Value','')
            app.tags.ShowInt='Single';
            app.tags.interpolated='n';
            app.tags.noised='n';
            app.tags.noisysim='n';
            app.tags.adjusted='n';
            app.tags.BaselineForInt='None';
            app.tags.ewrlsed='n';
            app.tags.ewrlsavged='n';
            app.tags.gvaled='n';
            app.dataexists='y';
            app.xlab='Field (mT)';
            set(app.UIAxes,'xdir','normal');
            set(app.gvaluesButton,'Value',0)
            set(app.FittingSliceLabel,"Visible","Off")
            
            set(app.StaticExpListBox,'Items',app.StaticExpStr)
            
            getexp(app)
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,1),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice))...
                max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Value changed function: Spinner
        function spinner(app, event)
            switch app.dataexists
                case 'n'
                    set(app.Spinner,'Value',1)
                    return
            end
            showselect=get(app.ShowButtonGroup,'SelectedObject');
            switch showselect.Text
                case 'Sim Only'
                    set(app.SimOnlyButton,'Value',0)
                    set(app.DataOnlyButton,'Value',0)
                    set(app.PlotBothButton,'Value',1)
            end
            app.currslice = app.Spinner.Value;
            sliceplts(app)
            
        end

        % Value changed function: gvaluesButton
        function gvalues(app, event)
            value = app.gvaluesButton.Value;
            showselect=get(app.ShowButtonGroup,'SelectedObject');
            selbutton=showselect.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            
            switch value
                case 1
                    explst=get(app.ExpListBox,'Items');
                    for i=1:numel(explst)
                        explst_f=strsplit(explst{i},'=');
                        switch explst_f{1}
                            case 'mwFreq'
                                mwFQ=str2double(explst_f{2});
                        end
                    end
                    if isnan(mwFQ)
                        set(app.gvaluesButton,'Value',0)
                        return
                    end
                    app.xlab='g values';
                    switch selbutton
                        case 'Sim Only'
                            switch app.simexists
                                case 'n'
                                    switch app.dataexists
                                        case 'y'
                                            gaxis=(((6.6261e-34)*mwFQ*1e9)./((9.2740e-24).*app.data.B)).*1e3;
                                            app.data.ogB=app.data.B;
                                            app.data.B=gaxis;
                                            set(app.SimOnlyButton,'Value',0)
                                            set(app.DataOnlyButton,'Value',1)
                                            set(app.PlotBothButton,'Value',0)
                                        case 'n'
                                            set(app.gvaluesButton,'Value',0)
                                            return
                                    end
                            end
                            simgaxis=(((6.6261e-34)*mwFQ*1e9)./((9.2740e-24).*app.data.simB)).*1e3;
                            app.data.ogsimB=app.data.simB;
                            app.data.simB=simgaxis;
                        otherwise
                            switch app.dataexists
                                case 'n'
                                    switch app.simexists
                                        case 'y'
                                            simgaxis=(((6.6261e-34)*mwFQ*1e9)./((9.2740e-24).*app.data.simB)).*1e3;
                                            app.data.ogsimB=app.data.simB;
                                            app.data.simB=simgaxis;
                                            set(app.SimOnlyButton,'Value',1)
                                            set(app.DataOnlyButton,'Value',0)
                                            set(app.PlotBothButton,'Value',0)
                                        case 'n'
                                            set(app.gvaluesButton,'Value',0)
                                            return
                                    end
                            end
                            gaxis=(((6.6261e-34)*mwFQ*1e9)./((9.2740e-24).*app.data.B)).*1e3;
                            app.data.ogB=app.data.B;
                            app.data.B=gaxis;
                    end
                    togglesysplts(app,onORoff,selbutton)
                    set(app.UIAxes,'xdir','reverse');
                    set(app.gvaluesButton,'Text','Field (mT)')
                    app.tags.gvaled='y';
                case 0
                    app.xlab='Field (mT)';
                    switch app.tags.gvaled
                        case 'n'
                            set(app.gvaluesButton,'Value',0)
                            return
                    end
                    switch selbutton
                        case 'Sim Only'
                            app.data.simB=app.data.ogsimB;
                        otherwise
                            app.data.B=app.data.ogB;
                    end
                    togglesysplts(app,onORoff,selbutton)
                    set(app.UIAxes,'xdir','normal');
                    set(app.gvaluesButton,'Text','g values')
                    app.tags.gvaled='n';
            end
            
        end

        % Value changed function: Sys1Switch
        function sys1selswitch(app, event)
            selectedButton = get(app.ShowButtonGroup,'SelectedObject');
            selbutton=selectedButton.Text;
            onORoff{1}=app.Sys1Switch.Value;
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
        end

        % Value changed function: Sys1ListBox
        function Sys1List(app, event)
            set(app.Sys1Edit,'Value',app.Sys1ListBox.Value)
            set(app.Vary1Edit,'Value',app.Sys1ListBox.Value)
        end

        % Value changed function: Sys1Edit
        function EditSys1(app, event)
            sysstr = app.Sys1Edit.Value;
            syslst=get(app.Sys1ListBox,'Items');
            set(app.Sys1Edit,'Value','')
            set(app.Vary1Edit,'Value','')
            if contains(sysstr,'import')
                app.inp.Sys1Str=syslst;
                sysstr_splt=strsplit(sysstr,' ');
                sysstrct=evalin('base',sysstr_splt{2});
                sys_params=fieldnames(sysstrct);
                [syslst_new] = getsysstructs(app,app.SysStr,sysstrct,sys_params);
            else
                switch sysstr
                    case 'reset'
                        app.inp.Sys1Str=syslst;
                        set(app.Sys1ListBox,'Items',app.SysStr)
                        set(app.Sys1Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Sys1ListBox,'Items',app.inp.Sys1Str)
                        return
                    case 'get Sys1'
                        return
                    case 'get Sys2'
                        app.inp.Sys1Str=syslst;
                        sys2copy=get(app.Sys2ListBox,'Items');
                        set(app.Sys1ListBox,'Items',sys2copy)
                        return
                    case 'get Sys3'
                        app.inp.Sys1Str=syslst;
                        sys2copy=get(app.Sys3ListBox,'Items');
                        set(app.Sys1ListBox,'Items',sys2copy)
                        return
                    case 'get Sys4'
                        app.inp.Sys1Str=syslst;
                        sys2copy=get(app.Sys4ListBox,'Items');
                        set(app.Sys1ListBox,'Items',sys2copy)
                        return
                    otherwise
                        [syslst_new]=sysadd(app,sysstr,syslst);
                        app.inp.Sys1Str=syslst_new;
                end
            end
            
            set(app.Sys1ListBox,'Items',syslst_new)
            
        end

        % Value changed function: Sys2Switch
        function sys2selswitch(app, event)
            selectedButton = get(app.ShowButtonGroup,'SelectedObject');
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=app.Sys2Switch.Value;
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
        end

        % Value changed function: Sys2ListBox
        function Sys2List(app, event)
            set(app.Sys2Edit,'Value',app.Sys2ListBox.Value)
            set(app.Vary2Edit,'Value',app.Sys2ListBox.Value)
        end

        % Value changed function: Sys2Edit
        function EditSys2(app, event)
            sysstr = app.Sys2Edit.Value;
            syslst=get(app.Sys2ListBox,'Items');
            set(app.Sys2Edit,'Value','')
            set(app.Vary2Edit,'Value','')
            if contains(sysstr,'import')
                app.inp.Sys2Str=syslst;
                sysstr_splt=strsplit(sysstr,' ');
                sysstrct=evalin('base',sysstr_splt{2});
                sys_params=fieldnames(sysstrct);
                [syslst_new] = getsysstructs(app,app.SysStr,sysstrct,sys_params);
            else
                switch sysstr
                    case 'reset'
                        app.inp.Sys2Str=syslst;
                        set(app.Sys2ListBox,'Items',app.SysStr)
                        set(app.Sys2Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Sys2ListBox,'Items',app.inp.Sys2Str)
                        return
                    case 'get Sys2'
                        return
                    case 'get Sys1'
                        app.inp.Sys2Str=syslst;
                        sys2copy=get(app.Sys1ListBox,'Items');
                        set(app.Sys2ListBox,'Items',sys2copy)
                        return
                    case 'get Sys3'
                        app.inp.Sys2Str=syslst;
                        sys2copy=get(app.Sys3ListBox,'Items');
                        set(app.Sys2ListBox,'Items',sys2copy)
                        return
                    case 'get Sys4'
                        app.inp.Sys2Str=syslst;
                        sys2copy=get(app.Sys4ListBox,'Items');
                        set(app.Sys2ListBox,'Items',sys2copy)
                        return
                    otherwise
                        [syslst_new]=sysadd(app,sysstr,syslst);
                        app.inp.Sys2Str=syslst_new;
                end
            end
            
            set(app.Sys2ListBox,'Items',syslst_new)
            
        end

        % Value changed function: Sys3Switch
        function sys3selswitch(app, event)
            selectedButton = get(app.ShowButtonGroup,'SelectedObject');
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=app.Sys3Switch.Value;
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
        end

        % Value changed function: Sys3ListBox
        function Sys3List(app, event)
            set(app.Sys3Edit,'Value',app.Sys3ListBox.Value)
            set(app.Vary3Edit,'Value',app.Sys3ListBox.Value)
        end

        % Value changed function: Sys3Edit
        function EditSys3(app, event)
            sysstr = app.Sys3Edit.Value;
            syslst=get(app.Sys3ListBox,'Items');
            set(app.Sys3Edit,'Value','')
            set(app.Vary3Edit,'Value','')
            if contains(sysstr,'import')
                app.inp.Sys3Str=syslst;
                sysstr_splt=strsplit(sysstr,' ');
                sysstrct=evalin('base',sysstr_splt{2});
                sys_params=fieldnames(sysstrct);
                [syslst_new] = getsysstructs(app,app.SysStr,sysstrct,sys_params);
            else
                switch sysstr
                    case 'reset'
                        app.inp.Sys3Str=syslst;
                        set(app.Sys3ListBox,'Items',app.SysStr)
                        set(app.Sys3Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Sys3ListBox,'Items',app.inp.Sys3Str)
                        return
                    case 'get Sys3'
                        return
                    case 'get Sys1'
                        app.inp.Sys3Str=syslst;
                        sys2copy=get(app.Sys1ListBox,'Items');
                        set(app.Sys3ListBox,'Items',sys2copy)
                        return
                    case 'get Sys2'
                        app.inp.Sys3Str=syslst;
                        sys2copy=get(app.Sys2ListBox,'Items');
                        set(app.Sys3ListBox,'Items',sys2copy)
                        return
                    case 'get Sys4'
                        app.inp.Sys3Str=syslst;
                        sys2copy=get(app.Sys4ListBox,'Items');
                        set(app.Sys3ListBox,'Items',sys2copy)
                        return
                    otherwise
                        [syslst_new]=sysadd(app,sysstr,syslst);
                        app.inp.Sys3Str=syslst_new;
                end
            end
            
            set(app.Sys3ListBox,'Items',syslst_new)
            
        end

        % Value changed function: Sys4Switch
        function sys4selswitch(app, event)
            selectedButton = get(app.ShowButtonGroup,'SelectedObject');
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=app.Sys4Switch.Value;
            togglesysplts(app,onORoff,selbutton)
        end

        % Value changed function: Sys4ListBox
        function Sys4List(app, event)
            set(app.Sys4Edit,'Value',app.Sys4ListBox.Value)
            set(app.Vary4Edit,'Value',app.Sys4ListBox.Value)
        end

        % Value changed function: Sys4Edit
        function EditSys4(app, event)
            sysstr = app.Sys4Edit.Value;
            syslst=get(app.Sys4ListBox,'Items');
            set(app.Sys4Edit,'Value','')
            set(app.Vary4Edit,'Value','')
            if contains(sysstr,'import')
                app.inp.Sys4Str=syslst;
                sysstr_splt=strsplit(sysstr,' ');
                sysstrct=evalin('base',sysstr_splt{2});
                sys_params=fieldnames(sysstrct);
                [syslst_new] = getsysstructs(app,app.SysStr,sysstrct,sys_params);
            else
                switch sysstr
                    case 'reset'
                        app.inp.Sys4Str=syslst;
                        set(app.Sys4ListBox,'Items',app.SysStr)
                        set(app.Sys4Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Sys4ListBox,'Items',app.inp.Sys4Str)
                        return
                    case 'get Sys4'
                        return
                    case 'get Sys1'
                        app.inp.Sys4Str=syslst;
                        sys2copy=get(app.Sys1ListBox,'Items');
                        set(app.Sys4ListBox,'Items',sys2copy)
                        return
                    case 'get Sys2'
                        app.inp.Sys4Str=syslst;
                        sys2copy=get(app.Sys2ListBox,'Items');
                        set(app.Sys4ListBox,'Items',sys2copy)
                        return
                    case 'get Sys3'
                        app.inp.Sys4Str=syslst;
                        sys2copy=get(app.Sys3ListBox,'Items');
                        set(app.Sys4ListBox,'Items',sys2copy)
                        return
                    otherwise
                        [syslst_new]=sysadd(app,sysstr,syslst);
                        app.inp.Sys4Str=syslst_new;
                end
            end
            
            set(app.Sys4ListBox,'Items',syslst_new)
            
        end

        % Value changed function: ExpListBox
        function ExpList(app, event)
            set(app.ExpEdit,'Value',app.ExpListBox.Value)
        end

        % Value changed function: ExpEdit
        function EditExp(app, event)
            parstr = app.ExpEdit.Value;
            parlst=get(app.ExpListBox,'Items');
            set(app.ExpEdit,'Value','')
            if contains(parstr,'import')
                app.inp.ExpStr=parlst;
                parstr_splt=strsplit(parstr,' ');
                expstrct=evalin('base',parstr_splt{2});
                exp_params=fieldnames(expstrct);
                [parlst_new] = getsysstructs(app,parlst,expstrct,exp_params);
            else
                switch parstr
                    case 'reset'
                        app.inp.ExpStr=parlst;
                        set(app.ExpListBox,'Items',app.ExpStr)
                        return
                    case 'revert'
                        set(app.ExpListBox,'Items',app.inp.ExpStr)
                        return
                    otherwise
                        [parlst_new]=paradd(app,parstr,parlst);
                        app.inp.ExpStr=parlst_new;
                end
            end
            
            set(app.ExpListBox,'Items',parlst_new)
            
        end

        % Button pushed function: FromDataButton
        function getexp(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            inexplst=get(app.ExpListBox,'Items');
            for i=1:numel(inexplst)
                explst_f=strsplit(inexplst{i},'=');
                switch explst_f{1}
                    case 'nPoints'
                        inexplst{i}=['nPoints','=',num2str(length(app.data.B))];
                    case 'mwFreq'
                        inexplst{i}=['mwFreq','=',app.params.MFQ];
                    case 'ModAmp'
                        inexplst{i}=['ModAmp','=',app.params.modA];
                    case 'Harmonic'
                        inexplst{i}=['Harmonic','=',app.params.Harm];
                    case 'Range'
                        inexplst{i}=['Range','=','[',num2str(min(app.data.B)),' ',num2str(max(app.data.B)),']'];
                    otherwise
                        inexplst{i}=inexplst{i};
                end
            end
            
            set(app.ExpListBox,'Items',inexplst)
            app.inp.Exp=inexplst;
            
            staticexplst=get(app.StaticExpListBox,'Items');
            for i=1:numel(staticexplst)
                explst_f=strsplit(staticexplst{i},'=');
                switch explst_f{1}
                    case 'Q Value'
                        staticexplst{i}=['Q Value','=',app.params.QVal];
                    case 'Time Constant'
                        staticexplst{i}=['Time Constant','=',app.params.TC, ' ms'];
                    case 'Conversion Time'
                        staticexplst{i}=['Conversion Time','=',app.params.CT,' ms'];
                    case 'Reciever Gain'
                        if str2double(app.params.RG)<61
                            staticexplst{i}=['Reciever Gain','=',app.params.RG,  ' dB'];
                        else
                            staticexplst{i}=['Reciever Gain','=',app.params.RG];
                        end
                    case 'Power (mW)'
                        staticexplst{i}=['Power','=',app.params.POW,' mW'];
                    case 'Power (dB)'
                        staticexplst{i}=['Power','=',app.params.ATTN,' dB'];
                    case 'Averages'
                        staticexplst{i}=['Averages','=',app.params.nScans];
                    case 'Temperature'
                        staticexplst{i}=['Temperature','=',app.params.Temp,' K'];
                    otherwise
                        staticexplst{i}=staticexplst{i};
                end
            end
            
            set(app.StaticExpListBox,'Items',staticexplst)
            
        end

        % Value changing function: FrequencySlider
        function freqChanging(app, event)
            
            switch app.simexists
                case 'y'
                    switch app.dataexists
                        case 'n'
                            return
                    end
                case 'n'
                    return
            end
            
            slide_freq(app,event.Value)
            
            set(app.FreqSpinner,'Value',event.Value)
            
        end

        % Value changed function: FreqSpinner, FrequencySlider
        function freqChanged(app, event)
            
            switch app.simexists
                case 'y'
                    switch app.dataexists
                        case 'n'
                            return
                    end
                case 'n'
                    return
            end
            
            slide_freq(app,event.Value)
            
            inexplst=get(app.ExpListBox,'Items');
            for i=1:numel(inexplst)
                explst_f=strsplit(inexplst{i},'=');
                switch explst_f{1}
                    case 'mwFreq'
                        inexplst{i}=['mwFreq','=',num2str(event.Value,'%2.6f')];
                    otherwise
                        inexplst{i}=inexplst{i};
                end
            end
            set(app.ExpListBox,'Items',inexplst)
            app.inp.ExpStr=inexplst;
            
        end

        % Value changed function: OptListBox
        function OptList(app, event)
            set(app.OptEdit,'Value',app.OptListBox.Value)
        end

        % Value changed function: OptEdit
        function EditOpt(app, event)
            parstr = app.OptEdit.Value;
            parlst=get(app.OptListBox,'Items');
            set(app.OptEdit,'Value','')
            if contains(parstr,'import')
                app.inp.OptStr=parlst;
                parstr_splt=strsplit(parstr,' ');
                optstrct=evalin('base',parstr_splt{2});
                opt_params=fieldnames(optstrct);
                [parlst_new] = getsysstructs(app,parlst,optstrct,opt_params);
            else
                switch parstr
                    case 'reset'
                        app.inp.OptStr=parlst;
                        set(app.OptListBox,'Items',app.OptStr)
                        return
                    case 'revert'
                        set(app.OptListBox,'Items',app.inp.OptStr)
                        return
                    otherwise
                        [parlst_new]=paradd(app,parstr,parlst);
                        app.inp.OptStr=parlst_new;
                end
            end
            
            set(app.OptListBox,'Items',parlst_new)
            
        end

        % Value changed function: Vary1ListBox
        function Vary1List(app, event)
            set(app.Vary1Edit,'Value',app.Vary1ListBox.Value)
        end

        % Value changed function: Vary1Edit
        function EditVary1(app, event)
            varystr = app.Vary1Edit.Value;
            varylst=get(app.Vary1ListBox,'Items');
            set(app.Vary1Edit,'Value','')
            if contains(varystr,'import')
                app.inp.Vary1Str=varylst;
                varystr_splt=strsplit(varystr,' ');
                varystrct=evalin('base',varystr_splt{2});
                vary_params=fieldnames(varystrct);
                for i=1:numel(vary_params)
                    vrylst{i}=[vary_params{i},'='];
                end
                [varylst_new] = getsysstructs(app,vrylst,varystrct,vary_params);
            else
                switch varystr
                    case 'reset'
                        app.inp.Vary1Str=varylst;
                        set(app.Vary1ListBox,'Items',{})
                        set(app.Vary1Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Vary1ListBox,'Items',app.inp.Vary1Str)
                        return
                    case 'get Vary1'
                        return
                    case 'get Vary2'
                        app.inp.Vary1Str=varylst;
                        vary2copy=get(app.Vary2ListBox,'Items');
                        set(app.Vary1ListBox,'Items',vary2copy)
                        return
                    case 'get Vary3'
                        app.inp.Vary1Str=varylst;
                        vary2copy=get(app.Vary3ListBox,'Items');
                        set(app.Vary1ListBox,'Items',vary2copy)
                        return
                    case 'get Vary4'
                        app.inp.Vary1Str=varylst;
                        vary2copy=get(app.Vary4ListBox,'Items');
                        set(app.Vary1ListBox,'Items',vary2copy)
                        return
                    otherwise
                        [varylst_new]=varyadd(app,varystr,varylst);
                        app.inp.Vary1Str=varylst_new;
                end
            end
            
            set(app.Vary1ListBox,'Items',varylst_new)
            set(app.Vary1Switch,'Value','On')
            
        end

        % Button pushed function: Vary1RemoveButton
        function Vary1Remove(app, event)
            varyitems=get(app.Vary1ListBox,'Items');
            varyval=get(app.Vary1ListBox,'Value');
            if numel(varyitems)==1
                new_varyitems={};
            else
                k=1;
                for i=1:numel(varyitems)-1
                    switch varyval
                        case varyitems{i}
                            t=i+1;
                            new_varyitems{i}=varyitems{t};
                            k=t+1;
                        otherwise
                            new_varyitems{i}=varyitems{k};
                            k=k+1;
                    end
                end
            end
            
            set(app.Vary1ListBox,'Items',new_varyitems)
            
        end

        % Value changed function: Vary2ListBox
        function Vary2List(app, event)
            set(app.Vary2Edit,'Value',app.Vary2ListBox.Value)
        end

        % Value changed function: Vary2Edit
        function EditVary2(app, event)
            varystr = app.Vary2Edit.Value;
            varylst=get(app.Vary2ListBox,'Items');
            set(app.Vary2Edit,'Value','')
            if contains(varystr,'import')
                app.inp.Vary2Str=varylst;
                varystr_splt=strsplit(varystr,' ');
                varystrct=evalin('base',varystr_splt{2});
                vary_params=fieldnames(varystrct);
                for i=1:numel(vary_params)
                    vrylst{i}=[vary_params{i},'='];
                end
                [varylst_new] = getsysstructs(app,vrylst,varystrct,vary_params);
            else
                switch varystr
                    case 'reset'
                        app.inp.Vary2Str=varylst;
                        set(app.Vary2ListBox,'Items',{})
                        set(app.Vary2Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Vary2ListBox,'Items',app.inp.Vary2Str)
                        return
                    case 'get Vary2'
                        return
                    case 'get Vary1'
                        app.inp.Vary2Str=varylst;
                        vary2copy=get(app.Vary1ListBox,'Items');
                        set(app.Vary2ListBox,'Items',vary2copy)
                        return
                    case 'get Vary3'
                        app.inp.Vary2Str=varylst;
                        vary2copy=get(app.Vary3ListBox,'Items');
                        set(app.Vary2ListBox,'Items',vary2copy)
                        return
                    case 'get Vary4'
                        app.inp.Vary2Str=varylst;
                        vary2copy=get(app.Vary4ListBox,'Items');
                        set(app.Vary2ListBox,'Items',vary2copy)
                        return
                    otherwise
                        [varylst_new]=varyadd(app,varystr,varylst);
                        app.inp.Vary2Str=varylst_new;
                end
            end
            
            set(app.Vary2ListBox,'Items',varylst_new)
            set(app.Vary2Switch,'Value','On')
            
        end

        % Button pushed function: Vary2RemoveButton
        function Vary2Remove(app, event)
            varyitems=get(app.Vary2ListBox,'Items');
            varyval=get(app.Vary2ListBox,'Value');
            if numel(varyitems)==1
                new_varyitems={};
            else
                k=1;
                for i=1:numel(varyitems)-1
                    switch varyval
                        case varyitems{i}
                            t=i+1;
                            new_varyitems{i}=varyitems{t};
                            k=t+1;
                        otherwise
                            new_varyitems{i}=varyitems{k};
                            k=k+1;
                    end
                end
            end
            
            set(app.Vary2ListBox,'Items',new_varyitems)
            
        end

        % Value changed function: Vary3ListBox
        function Vary3List(app, event)
            set(app.Vary3Edit,'Value',app.Vary3ListBox.Value)
        end

        % Value changed function: Vary3Edit
        function EditVary3(app, event)
            varystr = app.Vary3Edit.Value;
            varylst=get(app.Vary3ListBox,'Items');
            set(app.Vary3Edit,'Value','')
            if contains(varystr,'import')
                app.inp.Vary3Str=varylst;
                varystr_splt=strsplit(varystr,' ');
                varystrct=evalin('base',varystr_splt{2});
                vary_params=fieldnames(varystrct);
                for i=1:numel(vary_params)
                    vrylst{i}=[vary_params{i},'='];
                end
                [varylst_new] = getsysstructs(app,vrylst,varystrct,vary_params);
            else
                switch varystr
                    case 'reset'
                        app.inp.Vary3Str=varylst;
                        set(app.Vary3ListBox,'Items',{})
                        set(app.Vary3Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Vary3ListBox,'Items',app.inp.Vary3Str)
                        return
                    case 'get Vary3'
                        return
                    case 'get Vary2'
                        app.inp.Vary3Str=varylst;
                        vary2copy=get(app.Vary2ListBox,'Items');
                        set(app.Vary3ListBox,'Items',vary2copy)
                        return
                    case 'get Vary1'
                        app.inp.Vary3Str=varylst;
                        vary2copy=get(app.Vary1ListBox,'Items');
                        set(app.Vary3ListBox,'Items',vary2copy)
                        return
                    case 'get Vary4'
                        app.inp.Vary3Str=varylst;
                        vary2copy=get(app.Vary4ListBox,'Items');
                        set(app.Vary3ListBox,'Items',vary2copy)
                        return
                    otherwise
                        [varylst_new]=varyadd(app,varystr,varylst);
                        app.inp.Vary3Str=varylst_new;
                end
            end
            
            set(app.Vary3ListBox,'Items',varylst_new)
            set(app.Vary3Switch,'Value','On')
            
        end

        % Button pushed function: Vary3RemoveButton
        function Vary3Remove(app, event)
            varyitems=get(app.Vary3ListBox,'Items');
            varyval=get(app.Vary3ListBox,'Value');
            if numel(varyitems)==1
                new_varyitems={};
            else
                k=1;
                for i=1:numel(varyitems)-1
                    switch varyval
                        case varyitems{i}
                            t=i+1;
                            new_varyitems{i}=varyitems{t};
                            k=t+1;
                        otherwise
                            new_varyitems{i}=varyitems{k};
                            k=k+1;
                    end
                end
            end
            
            set(app.Vary3ListBox,'Items',new_varyitems)
            
        end

        % Value changed function: Vary4ListBox
        function Vary4List(app, event)
            set(app.Vary4Edit,'Value',app.Vary4ListBox.Value)
        end

        % Value changed function: Vary4Edit
        function EditVary4(app, event)
            varystr = app.Vary4Edit.Value;
            varylst=get(app.Vary4ListBox,'Items');
            set(app.Vary4Edit,'Value','')
            if contains(varystr,'import')
                app.inp.Vary4Str=varylst;
                varystr_splt=strsplit(varystr,' ');
                varystrct=evalin('base',varystr_splt{2});
                vary_params=fieldnames(varystrct);
                for i=1:numel(vary_params)
                    vrylst{i}=[vary_params{i},'='];
                end
                [varylst_new] = getsysstructs(app,vrylst,varystrct,vary_params);
            else
                switch varystr
                    case 'reset'
                        app.inp.Vary4Str=varylst;
                        set(app.Vary4ListBox,'Items',{})
                        set(app.Vary4Switch,'Value','Off')
                        return
                    case 'revert'
                        set(app.Vary4ListBox,'Items',app.inp.Vary4Str)
                        return
                    case 'get Vary4'
                        return
                    case 'get Vary2'
                        app.inp.Vary4Str=varylst;
                        vary2copy=get(app.Vary2ListBox,'Items');
                        set(app.Vary4ListBox,'Items',vary2copy)
                        return
                    case 'get Vary3'
                        app.inp.Vary4Str=varylst;
                        vary2copy=get(app.Vary3ListBox,'Items');
                        set(app.Vary4ListBox,'Items',vary2copy)
                        return
                    case 'get Vary1'
                        app.inp.Vary4Str=varylst;
                        vary2copy=get(app.Vary1ListBox,'Items');
                        set(app.Vary4ListBox,'Items',vary2copy)
                        return
                    otherwise
                        [varylst_new]=varyadd(app,varystr,varylst);
                        app.inp.Vary4Str=varylst_new;
                end
            end
            
            set(app.Vary4ListBox,'Items',varylst_new)
            set(app.Vary4Switch,'Value','On')
            
        end

        % Button pushed function: Vary4RemoveButton
        function Vary4Remove(app, event)
            varyitems=get(app.Vary4ListBox,'Items');
            varyval=get(app.Vary4ListBox,'Value');
            if numel(varyitems)==1
                new_varyitems={};
            else
                k=1;
                for i=1:numel(varyitems)-1
                    switch varyval
                        case varyitems{i}
                            t=i+1;
                            new_varyitems{i}=varyitems{t};
                            k=t+1;
                        otherwise
                            new_varyitems{i}=varyitems{k};
                            k=k+1;
                    end
                end
            end
            
            set(app.Vary4ListBox,'Items',new_varyitems)
            
        end

        % Value changed function: FitOptListBox
        function FitOptList(app, event)
            set(app.FitOptEdit,'Value',app.FitOptListBox.Value)
        end

        % Value changed function: FitOptEdit
        function EditFitOpt(app, event)
            parstr = app.FitOptEdit.Value;
            parlst=get(app.FitOptListBox,'Items');
            set(app.FitOptEdit,'Value','')
            if contains(parstr,'import')
                app.inp.FitOptStr=parlst;
                parstr_splt=strsplit(parstr,' ');
                fitoptstrct=evalin('base',parstr_splt{2});
                fitopt_params=fieldnames(fitoptstrct);
                [parlst_new] = getsysstructs(app,parlst,fitoptstrct,fitopt_params);
            else
                switch parstr
                    case 'reset'
                        app.inp.FitOptStr=parlst;
                        set(app.FitOptListBox,'Items',app.FitOptStr)
                        return
                    case 'revert'
                        set(app.FitOptListBox,'Items',app.inp.FitOptStr)
                        return
                    otherwise
                        [parlst_new]=paradd(app,parstr,parlst);
                        app.inp.FitOptStr=parlst_new;
                end
            end
            
            set(app.FitOptListBox,'Items',parlst_new)
            
        end

        % Button pushed function: LoadButton
        function LoadBackground(app, event)
            switch app.dataexists
                case 'n'
                    set(app.BackgroundFilenameLabel,'Text','Load data first')
                    return
            end
            app.cwEPRapp.Visible = 'Off';
            try
                [name,path]=uigetfile({'*.*';'*.DTA';'*.spc';'*.ESR';'*.txt';'*.dat';'*.xml'});
                [rootnm,exten]=strtok(name,'.');
                [~,userbase,~,~]=opendata(app,strcat(path,filesep,name),exten);
            catch
                app.cwEPRapp.Visible = 'On';
                return
            end
            app.cwEPRapp.Visible = 'On';
            
            if length(userbase)~=length(app.data.spc(:,1)) || length(userbase(1,:))>1
                set(app.BackgroundFilenameLabel,'Text','axes mismatch')
                return
            else
                app.data.userbase=[];
                app.data.userbase=basecorr(userbase,1,0);
                set(app.BackgroundFilenameLabel,'Text',rootnm)
            end
            
        end

        % Button pushed function: BackgroundButton
        function Background(app, event)
            try
                userbase=app.data.userbase;
            catch
                return
            end
            if ~isempty(app.data.userbase)
                app.data.subspc=[];
                app.data.fitln=[];
                for i=1:length(app.data.spc(1,:))
                    app.data.fitln(:,i)=userbase;
                    app.data.subspc(:,i)=app.data.spc(:,i)-userbase;
                end
                app.tags.baselined='y';
                app.tags.lastprocess='baseline';
                plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
                hold(app.UIAxes,'on')
                plot(app.UIAxes,app.data.B,app.data.subspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
                plot(app.UIAxes,app.data.B,app.data.fitln(:,app.currslice),'linewidth',app.plw,'color',app.pclr{2});
                xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
                ylim(app.UIAxes,[min(app.data.subspc(:,app.currslice)) max(app.data.subspc(:,app.currslice))].*1.05);
                ylabel(app.UIAxes,'Amplitude (a.u.)');
                xlabel(app.UIAxes,app.xlab);
                legend(app.UIAxes,'Original','Corrected','Background','location','northeast');
                hold(app.UIAxes,'off')
                set(app.DataOnlyButton,'Value',1)
                set(app.SimOnlyButton,'Value',0)
                set(app.PlotBothButton,'Value',0)
            end
            
        end

        % Callback function: PolyOrderSpinner, PolynomialFitButton,
        % SpectrumUsageSlider
        function PolynomialFit(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            app.data.subspc=[];
            app.data.fitln=[];
            polyorder=get(app.PolyOrderSpinner,'Value');
            baseperc=round(get(app.SpectrumUsageSlider,'Value'));
            l=length(app.data.spc(:,1));
            if baseperc == 0 || baseperc == 50
                l1 = 1;
                l2 = l;
            else
                l1=round(l*baseperc/100);
                l2=l-l1;
            end
            for i=1:length(app.data.spc(1,:))
                truncspc1=app.data.spc(1:l1,i);
                truncspc2=app.data.spc(l-l1+1:l,i);
                truncspc3=linspace(app.data.spc(l1,i),app.data.spc(l-l1+1,i),l-(2*l1))';
                truncspc(1:l1)=truncspc1;
                truncspc(l1+1:l-l1)=truncspc3;
                truncspc(l-l1+1:l)=truncspc2;
                [~,app.data.fitln(:,i)]=basecorr(truncspc',1,polyorder);
                app.data.subspc(:,i)=app.data.spc(:,i)-app.data.fitln(:,i);
            end
            app.tags.baselined='y';
            app.tags.lastprocess='baseline';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on')
            plot(app.UIAxes,app.data.B,app.data.subspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
            plot(app.UIAxes,app.data.B,app.data.fitln(:,app.currslice),'linewidth',app.plw,'color',app.pclr{2});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.subspc(:,app.currslice)) max(app.data.subspc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Original','Corrected','Baseline','location','northeast');
            if baseperc~=50
                plot(app.UIAxes,app.data.B(l1),0,'x','linewidth',4,'MarkerSize',16,'color','k')
                plot(app.UIAxes,app.data.B(l2),0,'x','linewidth',4,'MarkerSize',16,'color','k')
                legend(app.UIAxes,'Original','Corrected','Baseline',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northeast');
            end
            hold(app.UIAxes,'off')
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: msbackadjButton
        function msbackadjustment(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            app.data.subspc=[];
            app.data.fitln=[];
            warning Off
            %{
            inexplst=get(app.ExpListBox,'Items');
            for i=1:numel(inexplst)
                explst_f=strsplit(inexplst{i},'=');
                switch explst_f{1}
                    case 'nPoints'
                        npts=str2double(explst_f{2});
                    case 'ModAmp'
                        mdam=str2double(explst_f{2});
                end
            end
            if isempty(npts)
                npts=length(app.data.B);
            end
            if isempty(mdam)
                mdam=.1;
            end
            wf = @(B,npts,mdam)
            round((npts/(B(length(B))-B(1)))/(mdam*10));
            sf = @(B,npts,mdam) round(((npts/(B(length(B))-B(1)))/(mdam*10))/2);
            %}
            subspc = zeros(length(app.data.B),length(app.data.spc(1,:)));
            for i=1:length(app.data.spc(1,:))
                subspc(:,i) = msbackadj(app.data.B,app.data.spc(:,i));
                %                         'WindowSize',wf(app.data.B,npts,mdam),'StepSize',sf(app.data.B,npts,mdam),...
                %                         'RegressionMethod','pchip','EstimationMethod','quantile',...
                %                         'QuantileValue',0.1,'SmoothMethod','rloess','PreserveHeights',false);
                [app.data.subspc(:,i),~]=basecorr(subspc(:,i),1,0);
                app.data.fitln(:,i)=app.data.spc(:,i)-app.data.subspc(:,i);
            end
            warning On
            app.tags.baselined='y';
            app.tags.lastprocess='baseline';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on')
            plot(app.UIAxes,app.data.B,app.data.subspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
            plot(app.UIAxes,app.data.B,app.data.fitln(:,app.currslice),'linewidth',app.plw,'color',app.pclr{2});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.subspc(:,app.currslice)) max(app.data.subspc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Original','Baseline','Corrected','location','northeast');
            hold(app.UIAxes,'off')
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: ApplyButton
        function ApplyBaseline(app, event)
            switch app.tags.baselined
                case 'n'
                    return
            end
            app.data.org_subspc=app.data.spc;
            app.data.spc=[];
            switch get(app.NormalizeBaseCheckBox,'Value')
                case 1
                    app.data.spc=app.data.subspc./max(app.data.subspc(:));
                case 0
                    app.data.spc=app.data.subspc;
            end
            app.tags.subtracted='y';
            app.tags.baselined='n';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'New Spectrum','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: RevertBaseButton
        function RevertBaseline(app, event)
            switch app.tags.subtracted
                case 'y'
                    app.data.spc=[];
                    app.data.spc=app.data.org_subspc;
                case 'n'
                    switch app.tags.baselined
                        case 'n'
                            return
                    end
            end
            app.tags.subtracted='n';
            app.tags.baselined='n';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Callback function: ShowButton, StepsSpinner
        function ShowSmoothing(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            fsteps=get(app.StepsSpinner,'Value');
            if fsteps == 0
                app.RevertSmoothing()
                return
            end
            
            app.data.fspc=[];
            set(app.SmoothingBusyLabel,'Visible','On')
            drawnow
            
            selectedButton = get(app.SmoothingOptionsButtonGroup,'SelectedObject');
            selectedSmoothingButton=selectedButton.Text;
            switch selectedSmoothingButton
                case 'Moving Average'
                    app.data.fspc=datasmooth(app.data.spc,fsteps);
                case 'Binomial'
                    app.data.fspc=datasmooth(app.data.spc,fsteps,'binom');
                case 'Flat'
                    app.data.fspc=datasmooth(app.data.spc,fsteps,'flat');
                case 'Savitsky-Golay'
                    savgolstr=get(app.pdifEditField,'Value');
                    if isempty(savgolstr)
                        app.data.fspc=datasmooth(app.data.spc,fsteps,'savgol');
                    else
                        savgolvals=str2num(savgolstr);
                        if numel(savgolvals)==1
                            p=savgolvals(1);
                            app.data.fspc=datasmooth(app.data.spc,fsteps,'savgol',p);
                        elseif numel(savgolvals)==2
                            p=savgolvals(1);
                            dif=savgolvals(2);
                            app.data.fspc=datasmooth(app.data.spc,fsteps,'savgol',p,dif);
                        end
                    end
                case 'wdenoise'
                    wdlvl=str2double(get(app.wdOptionsEditField,'Value'));
                    if isempty(wdlvl) || isnan(wdlvl)
                        app.data.fspc=wdenoise(app.data.spc);
                    elseif wdlvl>floor(log2(length(app.data.B)))
                        wdlvl=floor(log2(length(app.data.B)));
                        set(app.wdOptionsEditField,'Value',num2str(wdlvl))
                        app.data.fspc=wdenoise(app.data.spc,wdlvl);
                    elseif ~isinteger(wdlvl) || wdlvl==0
                        wdlvl=round(wdlvl);
                        if wdlvl>floor(log2(length(app.data.B)))
                            wdlvl=floor(log2(length(app.data.B)));
                        elseif wdlvl==0
                            wdlvl=1;
                        end
                        set(app.wdOptionsEditField,'Value',num2str(wdlvl))
                        app.data.fspc=wdenoise(app.data.spc,wdlvl);
                    else
                        app.data.fspc=wdenoise(app.data.spc,wdlvl);
                    end
            end
            app.tags.filtered='y';
            app.tags.lastprocess='smoothing';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.data.B,app.data.fspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Original','Smoothed','location','northeast');
            hold(app.UIAxes,'off');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
            set(app.SmoothingBusyLabel,'Visible','Off')
            
        end

        % Button pushed function: SmoothingApplyButton
        function ApplySmoothing(app, event)
            switch app.tags.filtered
                case 'n'
                    return
            end
            app.data.org_fspc=[];
            app.data.org_fspc=app.data.spc;
            app.data.spc=[];
            app.data.spc=app.data.fspc;
            app.tags.smoothed='y';
            app.tags.lastprocess='none';
            app.tags.filtered='n';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'New Spectrum','location','northeast')
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: SmoothingRevertButton
        function RevertSmoothing(app, event)
            switch app.tags.smoothed
                case 'y'
                    app.data.spc=[];
                    app.data.spc=app.data.org_fspc;
                case 'n'
                    switch app.tags.filtered
                        case 'n'
                            return
                    end
            end
            app.tags.smoothed='n';
            app.tags.filtered='n';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            set(app.StepsSpinner, 'Value', 0)
            
        end

        % Callback function: BaselineCorrButtonGroup,
        % IntSpectrumUsageSlider, IntegrateButton
        function Integrate(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            
            baseperc=round(get(app.IntSpectrumUsageSlider,'Value'));
            
            [app.data.intgrted, app.data.dblintgrted, app.data.intgrl,...
                app.data.dblintgrl, l1, l2] = integrate_calculate(app, baseperc, app.data.spc);
            
            set(app.DIEditField,'Value',num2str(app.data.dblintgrl(app.currslice)))
            plt_fac_dbl = (max(app.data.spc(:,app.currslice))/max(app.data.dblintgrted(:,app.currslice)));
            plt_fac_sngl = (max(app.data.spc(:,app.currslice))/max(app.data.intgrted(:,app.currslice)));
            app.tags.integrated='y';
            app.tags.lastprocess='integral';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            showint=get(app.ShowIntButtonGroup,'SelectedObject');
            app.tags.ShowInt=showint.Text;
            switch app.tags.ShowInt
                case 'Single'
                    plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                    legend(app.UIAxes,'Data','First Integral','location','northeast');
                    plot(app.UIAxes,app.data.B(l1),app.data.spc(l1,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    plot(app.UIAxes,app.data.B(l2),app.data.spc(l2,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    legend(app.UIAxes,'Data','First Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                case 'Double'
                    plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{1});
                    legend(app.UIAxes,'Data','Double Integral','location','northwest');
                    plot(app.UIAxes,app.data.B(l1),app.data.spc(l1,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    plot(app.UIAxes,app.data.B(l2),app.data.spc(l2,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    legend(app.UIAxes,'Data','Double Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                case 'Both'
                    plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                    plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{2});
                    legend(app.UIAxes,'Original','First Integral','Double Integral','location','northwest');
                    plot(app.UIAxes,app.data.B(l1),app.data.spc(l1,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    plot(app.UIAxes,app.data.B(l2),app.data.spc(l2,app.currslice),'x','linewidth',4,'MarkerSize',16,'color','k')
                    legend(app.UIAxes,'Data','First Integral','Double Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                case 'Neither'
                    legend(app.UIAxes,'Data','location','northeast');
            end
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            hold(app.UIAxes,'off');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Value changed function: IntOrderSpinner
        function integrate_spinner(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            
            baseselect=get(app.BaselineCorrButtonGroup,'SelectedObject');
            switch baseselect.Text
                case "None"
                    return
                case "Polynomial, order:"
                    app.Integrate()
            end
            
        end

        % Button pushed function: AutoButton
        function auto_percent(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            set(app.AutoBusyLabel,'Visible','On')
            drawnow
            baseperc=2:1:49;
            dblintgrl = [];
            for j=1:length(app.data.spc(1,:))
                for k=1:length(baseperc)
                    [~,~,~,dblintgrl(j, k),~,~] = integrate_calculate(app, baseperc(k), app.data.spc(:,j));
                end
            end
            
            dblintgrl_sum = sum(dblintgrl, 1);
            [~,p,~]=find(dblintgrl_sum==max(dblintgrl_sum));
            set(app.IntSpectrumUsageSlider,'Value', baseperc(p))
            app.Integrate()
            set(app.AutoBusyLabel,'Visible','Off')
            
        end

        % Value changed function: TrapzCheckBox
        function trapz_select(app, event)
            switch app.TrapzCheckBox.Value
                case 1
                    set(app.SumCheckBox,'Value',0)
                case 0
                    set(app.SumCheckBox,'Value',1)
            end
            app.Integrate(app)
            
        end

        % Value changed function: SumCheckBox
        function sum_select(app, event)
            switch app.SumCheckBox.Value
                case 1
                    set(app.TrapzCheckBox,'Value',0)
                case 0
                    set(app.TrapzCheckBox,'Value',1)
            end
            app.Integrate(app)
            
        end

        % Selection changed function: ShowIntButtonGroup
        function ShowIntegrals(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            selectedButton = app.ShowIntButtonGroup.SelectedObject;
            app.tags.ShowInt=selectedButton.Text;
            baseperc = round(app.IntSpectrumUsageSlider.Value);
            spc=app.data.spc(:,app.currslice);
            l=length(app.data.B);
            if baseperc == 0 || baseperc == 50
                l1=1;
                l2=l;
            else
                l1=round(l*baseperc/100);
                l2=l-l1;
            end
            plt_fac_dbl = (max(spc)/max(app.data.dblintgrted(:,app.currslice)));
            plt_fac_sngl = (max(spc)/max(app.data.intgrted(:,app.currslice)));
            switch app.tags.integrated
                case 'y'
                    plot(app.UIAxes,app.data.B,spc,'linewidth',app.dlw,'color',app.dclr);
                    hold(app.UIAxes,'on');
                    switch app.tags.ShowInt
                        case 'Single'
                            plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                            legend(app.UIAxes,'Data','First Integral','location','northeast');
                            plot(app.UIAxes,app.data.B(l1),spc(l1),'x','linewidth',4,'MarkerSize',16,'color','k')
                            plot(app.UIAxes,app.data.B(l2),spc(l2),'x','linewidth',4,'MarkerSize',16,'color','k')
                            legend(app.UIAxes,'Data','First Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                        case 'Double'
                            plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{1});
                            legend(app.UIAxes,'Data','Double Integral','location','northwest');
                            plot(app.UIAxes,app.data.B(l1),spc(l1),'x','linewidth',4,'MarkerSize',16,'color','k')
                            plot(app.UIAxes,app.data.B(l2),spc(l2),'x','linewidth',4,'MarkerSize',16,'color','k')
                            legend(app.UIAxes,'Data','Double Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                        case 'Both'
                            plot(app.UIAxes,app.data.B,app.data.intgrted(:,app.currslice).*plt_fac_sngl,'linewidth',app.plw,'color',app.pclr{1});
                            plot(app.UIAxes,app.data.B,app.data.dblintgrted(:,app.currslice).*plt_fac_dbl,'linewidth',app.plw,'color',app.pclr{2});
                            legend(app.UIAxes,'Original','First Integral','Double Integral','location','northwest');
                            plot(app.UIAxes,app.data.B(l1),spc(l1),'x','linewidth',4,'MarkerSize',16,'color','k')
                            plot(app.UIAxes,app.data.B(l2),spc(l2),'x','linewidth',4,'MarkerSize',16,'color','k')
                            legend(app.UIAxes,'Data','First Integral','Double Integral',[num2str(round(app.data.B(l1),1)),' mT'],[num2str(round(app.data.B(l2),1)),' mT'],'location','northwest');
                        case 'Neither'
                            legend(app.UIAxes,'Data','location','northeast');
                    end
                    ylim(app.UIAxes,[min(spc) max(spc)].*1.05);
                    xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
                    ylabel(app.UIAxes,'Amplitude (a.u.)');
                    xlabel(app.UIAxes,app.xlab);
                    hold(app.UIAxes,'off');
                    set(app.DataOnlyButton,'Value',1)
                    set(app.SimOnlyButton,'Value',0)
                    set(app.PlotBothButton,'Value',0)
            end
        end

        % Callback function: InterpolateButton, PointsSpinner
        function Interpolate(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            oldlngth=length(app.data.B);
            newlngth=get(app.PointsSpinner,'Value');
            if isempty(newlngth) || newlngth==oldlngth || newlngth==0
                set(app.PointsSpinner,'Value',oldlngth)
                return
            end
            app.data.interpB=[];
            app.data.interpspc=[];
            newB=(linspace(app.data.B(1),app.data.B(oldlngth),newlngth))';
            app.data.interpB=newB;
            for j=1:length(app.data.spc(1,:))
                app.data.interpspc(:,j)=interp1(app.data.B,app.data.spc(:,j),newB,'pchip','extrap');
            end
            set(app.PointsSpinner,'Value',newlngth)
            app.tags.adjusted='interp';
            app.tags.lastprocess='interpolated';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.data.interpB,app.data.interpspc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
            legend(app.UIAxes,'Original','Interpolated','location','northeast');
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            hold(app.UIAxes,'off');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Callback function: AddNoiseButton, NoiseModelButtonGroup,
        % SNRSpinner
        function AddNoise(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            SNR=get(app.SNRSpinner,'Value');
            if isempty(SNR)
                return
            end
            app.data.noisespc=[];
            model=get(app.NoiseModelButtonGroup,'SelectedObject');
            switch model.Text
                case '1/f'
                    mdl='f';
                case 'Uniform'
                    mdl='u';
                case 'Gaussian'
                    mdl='n';
            end
            app.data.noisespc=addnoise(app.data.spc,SNR,mdl);
            app.tags.adjusted='noised';
            app.tags.lastprocess='noised';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.data.B,app.data.noisespc(:,app.currslice),'linewidth',app.plw,'color',app.pclr{1});
            legend(app.UIAxes,'Original','Noised','location','northeast');
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.noisespc(:,app.currslice)) max(app.data.noisespc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            hold(app.UIAxes,'off');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: NoisySimButton
        function NoisySim(app, event)
            switch app.simexists
                case 'n'
                    return
            end
            SNR=get(app.SNRSpinner,'Value');
            if isempty(SNR)
                return
            end
            
            show=get(app.ShowButtonGroup,'SelectedObject');
            switch show.Text
                case 'Sim Only'
                    model=get(app.NoiseModelButtonGroup,'SelectedObject');
                    switch model.Text
                        case '1/f'
                            mdl='f';
                        case 'Uniform'
                            mdl='u';
                        case 'Gaussian'
                            mdl='n';
                    end
                    simspc=sum(app.data.simspc,2);
                    noisespc=addnoise(simspc,SNR,mdl);
                    try
                        noisesimspc=rescaledata(noisespc,'maxabs');
                    catch
                        noisesimspc=rescale(noisespc,'maxabs');
                    end
                    plot(app.UIAxes,app.data.simB,noisesimspc,'linewidth',app.slw,'color',app.sclr{5});
                    xlim(app.UIAxes,[min(app.data.simB) max(app.data.simB)]);
                    ylim(app.UIAxes,[min(noisesimspc) max(noisesimspc)].*1.05);
                    ylabel(app.UIAxes,'Amplitude (a.u.)');
                    xlabel(app.UIAxes,app.xlab);
                    switch get(app.NoiseDTADSCCheckBox,'Value')
                        case 1
                            currdir=cd;
                            app.cwEPRapp.Visible = 'Off';
                            try
                                [fname,pth] = uiputfile({'*.*'});
                            catch
                                app.cwEPRapp.Visible = 'On';
                                return
                            end
                            app.cwEPRapp.Visible = 'On';
                            cd(pth)
                            eprsave(fname,app.data.simB,noisesimspc,fname,app.params.simFq);
                            cd(currdir)
                    end
                    switch get(app.NoiseasciiCheckBox,'Value')
                        case 1
                            svdata=zeros(length(app.data.simB),2);
                            svdata(:,1)=app.data.simB;
                            svdata(:,2:end)=noisesimspc;
                            currdir=cd;
                            app.cwEPRapp.Visible = 'Off';
                            try
                                [fnm,pth] = uiputfile({'*'});
                            catch
                                app.cwEPRapp.Visible = 'On';
                                return
                            end
                            app.cwEPRapp.Visible = 'On';
                            cd(pth)
                            fname=[fnm,'.dat'];
                            writematrix(svdata,fname);
                            cd(currdir)
                    end
                    switch get(app.NoisetoworkspaceCheckBox,'Value')
                        case 1
                            assignin('base','simB',app.data.simB)
                            assignin('base','noisySim',noisespc)
                    end
                    switch get(app.NoisetodataCheckBox,'Value')
                        case 1
                            app.currslice=1;
                            set(app.Spinner,'Visible','Off','Value',1)
                            set(app.SliceLabel,'Visible','Off')
                            app.tags.lastprocess='none';
                            app.tags.baselined='n';
                            app.tags.subtracted='n';
                            app.tags.filtered='n';
                            app.tags.smoothed='n';
                            app.tags.integrated='n';
                            app.tags.ShowInt='Single';
                            app.tags.interpolated='n';
                            app.tags.noised='n';
                            app.tags.adjusted='n';
                            app.tags.BaselineForInt='None';
                            app.tags.ewrlsed='n';
                            app.tags.ewrlsavged='n';
                            app.tags.secondmomented='n';
                            app.xlab='Field (mT)';
                            set(app.UIAxes,'xdir','normal');
                            set(app.gvaluesButton,'Value',0)
                            app.data.B=app.data.simB;
                            app.data.spc=noisespc;
                            switch app.dataexists
                                case 'n'
                                    app.data.orgB=app.data.B;
                                    app.data.orgspc=app.data.spc;
                                    app.dataexists='y';
                            end
                            legend(app.UIAxes,'New Data','location','northeast');
                        otherwise
                            legend(app.UIAxes,'Noisy Simulation','location','northeast');
                            app.tags.noisysim='y';
                    end
            end
            
        end

        % Value changed function: NoiseDTADSCCheckBox
        function NoiseDTADSC(app, event)
            switch app.NoiseDTADSCCheckBox.Value
                case 1
                    set(app.NoiseasciiCheckBox,'Value',0)
                    set(app.NoisetodataCheckBox,'Value',0)
            end
            
        end

        % Value changed function: NoiseasciiCheckBox
        function Noisetoascii(app, event)
            switch app.NoiseasciiCheckBox.Value
                case 1
                    set(app.NoiseDTADSCCheckBox,'Value',0)
                    set(app.NoisetodataCheckBox,'Value',0)
            end
        end

        % Value changed function: NoisetodataCheckBox
        function Noisetodata(app, event)
            switch app.NoisetodataCheckBox.Value
                case 1
                    set(app.NoiseasciiCheckBox,'Value',0)
                    set(app.NoiseDTADSCCheckBox,'Value',0)
                    set(app.NoisetoworkspaceCheckBox,'Value',0)
            end
            
        end

        % Button pushed function: ApplyAdjustmentsButton
        function ApplyAdjustments(app, event)
            switch app.tags.adjusted
                case 'interp'
                case 'noised'
                otherwise
                    return
            end
            app.data.org_aspc=app.data.spc;
            app.data.org_aB=app.data.B;
            switch app.tags.adjusted
                case 'interp'
                    app.tags.interpolated='y';
                    app.data.spc=[];
                    app.data.B=[];
                    app.data.spc=app.data.interpspc;
                    app.data.B=app.data.interpB;
                    newlngth=length(app.data.interpB);
                    explst=get(app.ExpListBox,'Items');
                    inexplst=explst;
                    for i=1:numel(explst)
                        explst_f=strsplit(explst{i},'=');
                        switch explst_f{1}
                            case 'nPoints'
                                inexplst{i}=['nPoints','=',num2str(newlngth)];
                            otherwise
                                inexplst{i}=inexplst{i};
                        end
                    end
                    set(app.ExpListBox,'Items',inexplst)
                case 'noised'
                    app.tags.noised='y';
                    app.data.spc=[];
                    app.data.spc=app.data.noisespc;
            end
            app.tags.adjusted='y';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: RevertAdjustmentsButton
        function RevertAdjustments(app, event)
            switch app.tags.noisysim
                case 'y'
                    onORoff{1}=get(app.Sys1Switch,'Value');
                    onORoff{2}=get(app.Sys2Switch,'Value');
                    onORoff{3}=get(app.Sys3Switch,'Value');
                    onORoff{4}=get(app.Sys4Switch,'Value');
                    showselect=get(app.ShowButtonGroup,'SelectedObject');
                    selbutton=showselect.Text;
                    showplots(app,selbutton,app.whchS,onORoff)
                    return
            end
            switch app.tags.adjusted
                case 'y'
                    app.data.spc=[];
                    app.data.B=[];
                    app.data.spc=app.data.org_aspc;
                    app.data.B=app.data.org_aB;
                case 'n'
                    switch app.tags.interpolated
                        case 'n'
                            switch app.tags.noised
                                case 'n'
                                    return
                            end
                    end
            end
            orglngth=length(app.data.orgB);
            set(app.PointsSpinner,'Value',orglngth)
            explst=get(app.ExpListBox,'Items');
            inexplst=explst;
            for i=1:numel(explst)
                explst_f=strsplit(explst{i},'=');
                switch explst_f{1}
                    case 'nPoints'
                        inexplst{i}=['nPoints','=',num2str(orglngth)];
                    otherwise
                        inexplst{i}=inexplst{i};
                end
            end
            set(app.ExpListBox,'Items',inexplst)
            app.tags.adjusted='n';
            app.tags.interpolated='n';
            app.tags.noised='n';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: ResetDataButton
        function ResetProcessing(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            app.data.spc=[];
            app.data.spc=app.data.orgspc;
            if length(app.data.spc(1,:))>1
                set(app.Spinner,'Visible','On')
            end
            app.data.B=[];
            app.data.B=app.data.orgB;
            orglngth=length(app.data.orgB);
            set(app.PointsSpinner,'Value',orglngth)
            explst=get(app.ExpListBox,'Items');
            inexplst=explst;
            for i=1:numel(explst)
                explst_f=strsplit(explst{i},'=');
                switch explst_f{1}
                    case 'nPoints'
                        inexplst{i}=['nPoints','=',num2str(orglngth)];
                    otherwise
                        inexplst{i}=inexplst{i};
                end
            end
            set(app.ExpListBox,'Items',inexplst)
            app.tags.interpolated='n';
            app.tags.noised='n';
            app.tags.adjusted='n';
            app.tags.filtered='n';
            app.tags.smoothed='n';
            app.tags.baselined='n';
            app.tags.subtracted='n';
            app.tags.lastprocess='none';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
        end

        % Button pushed function: ConversionToolButton
        function Convert(app, event)
            eprconvert
        end

        % Button pushed function: PeriodicTableButton
        function PeriodicTable(app, event)
            isotopes
        end

        % Value changed function: ConstantsDropDown
        function Constants(app, event)
            
            switch app.ConstantsDropDown.Value
                case 'Constants...'
                    set(app.ConstantLabel,'Text','')
                case 'amu'
                    set(app.ConstantLabel,'Text','1.6605e-27 kg')
                case 'angstrom'
                    set(app.ConstantLabel,'Text','1.0e-10 m')
                case 'avogadro'
                    set(app.ConstantLabel,'Text','6.0221e+23 /mol')
                case 'barn'
                    set(app.ConstantLabel,'Text','1.0e-28 m^2')
                case 'bmagn'
                    set(app.ConstantLabel,'Text','9.2740e-24 J/T')
                case 'bohrrad'
                    set(app.ConstantLabel,'Text','5.2918e-11 m')
                case 'boltzm'
                    set(app.ConstantLabel,'Text','1.3806e-23 J/K')
                case 'clight'
                    set(app.ConstantLabel,'Text','299792458 m/s')
                case 'echarge'
                    set(app.ConstantLabel,'Text','1.6022e-19 C')
                case 'emass'
                    set(app.ConstantLabel,'Text','9.1094e-31 kg')
                case 'eps0'
                    set(app.ConstantLabel,'Text','8.8542e-12 C/V*M')
                case 'evolt'
                    set(app.ConstantLabel,'Text','1.6022e-19 J')
                case 'faraday'
                    set(app.ConstantLabel,'Text','9.6485e+04 C/mol')
                case 'gfree'
                    set(app.ConstantLabel,'Text','2.0023')
                case 'hartree'
                    set(app.ConstantLabel,'Text','4.3597e-18 J')
                case 'hbar'
                    set(app.ConstantLabel,'Text','1.0546e-34 J*s')
                case 'molgas'
                    set(app.ConstantLabel,'Text','8.3145 J/K*mol')
                case 'mu0'
                    set(app.ConstantLabel,'Text','1.2566e-06 N/A^2')
                case 'nmagn'
                    set(app.ConstantLabel,'Text','5.0508e-27 J/T')
                case 'nmass'
                    set(app.ConstantLabel,'Text','1.6749e-27 kg')
                case 'planck'
                    set(app.ConstantLabel,'Text','6.6261e-34 J*s')
                case 'pmass'
                    set(app.ConstantLabel,'Text','1.6726e-27 kg')
                case 'rydberg'
                    set(app.ConstantLabel,'Text','1.0974e+07 /m')
                case 'pi'
                    set(app.ConstantLabel,'Text','3.14159265359')
            end
            
        end

        % Button pushed function: ndMomentButton
        function secondMoment(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            spc=app.data.spc;
            app.data.second_moment=[];
            polyorder=get(app.IntOrderSpinner,'Value');
            if isempty(polyorder)
                polyorder=0;
            end
            baseperc=round(get(app.IntSpectrumUsageSlider,'Value'));
            l=length(app.data.B(:,1));
            if baseperc == 0 || baseperc == 50
                l1 = 1;
                l2 = l;
            else
                l1=round(l*baseperc/100);
                l2=l-l1;
            end
            mB=mean(app.data.B);
            b_inc = mean(diff(app.data.B)*10);
            baseselect=get(app.BaselineCorrButtonGroup,'SelectedObject');
            app.data.second_moment=zeros(1,length(app.data.spc(1,:)));
            for i=1:length(spc(1,:))
                switch app.TrapzCheckBox.Value
                    case 1
                        frstint=cumtrapz(spc(:,i)).*b_inc;
                    case 0
                        frstint=cumsum(spc(:,i)).*b_inc;
                end
                truncspc1=frstint(1:l1);
                truncspc2=frstint(l-l1+1:l);
                truncspc3=linspace(frstint(l1),frstint(l-l1+1),l-(2*l1))';
                truncspc(1:l1)=truncspc1;
                truncspc(l1+1:l-l1)=truncspc3;
                truncspc(l-l1+1:l)=truncspc2;
                switch baseselect.Text
                    case 'None'
                        spc2=frstint;
                    case 'Polynomial, order:'
                        [~,fitln]=basecorr(truncspc',1,polyorder);
                        spc2=frstint-fitln;
                end
                secmom=((app.data.B-mB).^2).*spc2;
                switch app.TrapzCheckBox.Value
                    case 1
                        app.data.second_moment(i)=round(trapz(secmom(l1:l2))*b_inc);
                    case 0
                        app.data.second_moment(i)=round(sum(secmom(l1:l2))*b_inc);
                end
                %app.data.second_moment(i)=round(sqrt(second_m/length(second_m)));
            end
            app.tags.secondmomented='y';
            set(app.SecondMomentEdit,'Value',num2str(app.data.second_moment(:,app.currslice)))
            
        end

        % Button pushed function: ewrlsButton
        function ewrls(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            if length(app.data.spc(1,:))==1
                return
            end
            app.data.ewrls_spc=[];
            g_spl_=get(app.ewrlsoptions,'Value');
            set(app.ewrlsBusyLabel,'Visible','On')
            drawnow
            if contains(g_spl_,' ')
                g_spl=strsplit(g_spl_,' ');
            elseif contains(g_spl_,',')
                g_spl=strsplit(g_spl_,',');
            else
                g_spl=g_spl_;
            end
            if isempty(g_spl)
                p=30;
                lambda=0.96;
                nPreAv=0;
                delt=100;
                dir='fb';
            elseif ~iscell(g_spl)
                p=str2double(g_spl);
                lambda=0.96;
                nPreAv=0;
                delt=100;
                dir='fb';
            end
            if iscell(g_spl) && numel(g_spl)==2
                p=str2double(g_spl{1});
                lambda=str2double(g_spl{2});
                nPreAv=0;
                delt=100;
                dir='fb';
            elseif iscell(g_spl) && numel(g_spl)==3
                p=str2double(g_spl{1});
                lambda=str2double(g_spl{2});
                nPreAv=str2double(g_spl{3});
                delt=100;
                dir='fb';
            elseif iscell(g_spl) && numel(g_spl)==4
                p=str2double(g_spl{1});
                lambda=str2double(g_spl{2});
                nPreAv=str2double(g_spl{3});
                delt=str2double(g_spl{4});
                dir='fb';
            elseif iscell(g_spl) && numel(g_spl)==5
                p=str2double(g_spl{1});
                lambda=str2double(g_spl{2});
                nPreAv=str2double(g_spl{3});
                delt=str2double(g_spl{4});
                dir=g_spl{5};
            end
            app.data.ewrls_spc=ewrls(app.data.spc,p,lambda,nPreAv,delt,dir);
            app.tags.ewrlsed='y';
            set(app.Spinner,'Visible','Off')
            set(app.SliceLabel,'Visible','Off')
            set(app.ewrlsSpinner,'Limits',[1 length(app.data.spc(1,:))],'Value',1)
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.data.B,app.data.ewrls_spc,'linewidth',app.plw,'color',app.pclr{1});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.ewrls_spc) max(app.data.ewrls_spc)].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Averaged Spectrum','location','northeast');
            hold(app.UIAxes,'off');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
            set(app.ewrlsBusyLabel,'Visible','Off')
            
        end

        % Value changed function: ewrlsSpinner
        function ewrlsCompare(app, event)
            switch app.tags.ewrlsed
                case 'n'
                    set(app.ewrlsSpinner,'Value',1)
                    return
            end
            currentslice = app.ewrlsSpinner.Value;
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,currentslice),'linewidth',app.dlw,'color',app.dclr);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.data.B,app.data.ewrls_spc,'linewidth',app.plw,'color',app.pclr{1});
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.ewrls_spc) max(app.data.ewrls_spc)].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Original','ewrls','location','northeast');
            hold(app.UIAxes,'off');
            
        end

        % Button pushed function: ApplyewrlsButton
        function Applyewrls(app, event)
            switch app.tags.ewrlsed
                case 'n'
                    return
            end
            app.data.org_ewrlsspc=app.data.spc;
            app.data.spc=app.data.ewrls_spc;
            app.currslice=1;
            set(app.Spinner,'Value',1)
            app.tags.ewrlsavged='y';
            app.tags.ewrlsed='n';
            
            plot(app.UIAxes,app.data.B,app.data.spc,'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc) max(app.data.spc)].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'New 1D Spectrum','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
            set(app.ewrlsedLabel,'Visible','On')
            
        end

        % Button pushed function: RevertewrlsButton
        function Revertewrls(app, event)
            switch app.tags.ewrlsavged
                case 'y'
                    app.data.spc=app.data.org_ewrlsspc;
                    set(app.Spinner,'Visible','On')
                    set(app.SliceLabel,'Visible','On')
                case 'n'
                    switch app.tags.ewrlsed
                        case 'y'
                            set(app.Spinner,'Visible','On')
                            set(app.SliceLabel,'Visible','On')
                        case 'n'
                            return
                    end
            end
            
            app.tags.ewrlsavged='n';
            app.tags.ewrlsed='n';
            
            plot(app.UIAxes,app.data.B,app.data.spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
            xlim(app.UIAxes,[min(app.data.B) max(app.data.B)]);
            ylim(app.UIAxes,[min(app.data.spc(:,app.currslice)) max(app.data.spc(:,app.currslice))].*1.05);
            ylabel(app.UIAxes,'Amplitude (a.u.)');
            xlabel(app.UIAxes,app.xlab);
            legend(app.UIAxes,'Data','location','northeast');
            set(app.DataOnlyButton,'Value',1)
            set(app.SimOnlyButton,'Value',0)
            set(app.PlotBothButton,'Value',0)
            
            set(app.ewrlsedLabel,'Visible','Off')
            set(app.ewrlsBusyLabel,'Visible','Off')
            
        end

        % Value changed function: PlayButton
        function Play2D(app, event)
            switch app.dataexists
                case 'n'
                    set(app.PlayButton,'Value',0)
                    return
            end
            if length(app.data.spc(1,:))==1
                set(app.PlayButton,'Value',0)
                return
            end
            switch app.tags.ewrlsed
                case 'y'
                    set(app.PlayButton,'Value',0)
                    return
            end
            set(app.PlayButton,'Text','Stop')
            pauseval=get(app.PausesEditField,'Value');
            for i=1:length(app.data.spc(1,:))
                plystp=get(app.PlayButton,'Value');
                if plystp==0
                    sliceplts(app)
                    set(app.Spinner,'Value',app.currslice)
                    set(app.PlayButton,'Text','Play','Value',0)
                    break
                else
                    set(app.Spinner,'Value',i)
                    app.currslice=i;
                    sliceplts(app)
                    pause(pauseval)
                end
                hold(app.UIAxes,'off');
            end
            set(app.PlayButton,'Text','Play','Value',0)
            
        end

        % Button pushed function: SimulateButton
        function Simulate(app, event)
            tic
            switch app.xlab
                case 'g values'
                    app.data.B=app.data.ogB;
                    app.xlab='Field (mT)';
                    set(app.UIAxes,'xdir','normal');
                    set(app.gvaluesButton,'Value',0)
            end
            
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            funcval=get(app.FunctionDropDown,'Value');
            methodval=get(app.SimMethodDropDown,'Value');
            
            explst=get(app.ExpListBox,'Items');
            Exp=getexps(app,explst);
            
            if ~isfield(Exp,'mwFreq')
                set(app.SimTimeLabel,'Text','give mwFreq')
                return
            end
            app.params.simFq=Exp.('mwFreq');
            
            fq = app.params.simFq;
            set(app.FreqSpinner,'Value',fq)
            set(app.FrequencySlider,'Limits',[fq*.995 fq*1.005],'Value',fq)
            
            optlst=get(app.OptListBox,'Items');
            Opt = getopts(app,optlst);
            
            Sys={};
            s=0;
            switch onORoff{1}
                case 'On'
                    syslst1=get(app.Sys1ListBox,'Items');
                    s=s+1;
                    [Sys1,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst1);
                    if isempty(Sys1)
                        set(app.Sys1Switch,'Value','Off')
                        Sys1=[];
                        includeS(1)=0;
                    else
                        includeS(1)=1;
                    end
                case 'Off'
                    Sys1=[];
                    includeS(1)=0;
            end
            
            switch onORoff{2}
                case 'On'
                    syslst2=get(app.Sys2ListBox,'Items');
                    s=s+1;
                    [Sys2,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst2);
                    if isempty(Sys2)
                        set(app.Sys2Switch,'Value','Off')
                        Sys2=[];
                        includeS(2)=0;
                    else
                        includeS(2)=1;
                    end
                case 'Off'
                    Sys2=[];
                    includeS(2)=0;
            end
            
            switch onORoff{3}
                case 'On'
                    syslst3=get(app.Sys3ListBox,'Items');
                    s=s+1;
                    [Sys3,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst3);
                    if isempty(Sys3)
                        set(app.Sys3Switch,'Value','Off')
                        Sys3=[];
                        includeS(3)=0;
                    else
                        includeS(3)=1;
                    end
                case 'Off'
                    Sys3=[];
                    includeS(3)=0;
            end
            
            switch onORoff{4}
                case 'On'
                    syslst4=get(app.Sys4ListBox,'Items');
                    s=s+1;
                    [Sys4,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst4);
                    if isempty(Sys4)
                        set(app.Sys4Switch,'Value','Off')
                        Sys4=[];
                        includeS(4)=0;
                    else
                        includeS(4)=1;
                    end
                case 'Off'
                    Sys4=[];
                    includeS(4)=0;
            end
            if sum(includeS)==0
                set(app.FittingSliceLabel,'Visible','On','Text','No spin systems')
                return
            else
                set(app.FittingSliceLabel,'Visible','Off')
            end
            set(app.SimTimeLabel,'Text','busy...         ')
            drawnow
            
            if includeS(1)==1 && includeS(2)==0 && includeS(3)==0 && includeS(4)==0
                Sys=Sys1;
                app.whchS=1;
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==0 && includeS(4)==0
                Sys=Sys2;
                app.whchS=2;
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==1 && includeS(4)==0
                Sys=Sys3;
                app.whchS=3;
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==0 && includeS(4)==1
                Sys=Sys4;
                app.whchS=4;
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==0 && includeS(4)==0
                Sys={Sys1,Sys2};
                app.whchS=[1 2];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==1 && includeS(4)==0
                Sys={Sys1,Sys3};
                app.whchS=[1 3];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==0 && includeS(4)==1
                Sys={Sys1,Sys4};
                app.whchS=[1 4];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==1 && includeS(4)==0
                Sys={Sys2,Sys3};
                app.whchS=[2 3];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==0 && includeS(4)==1
                Sys={Sys2,Sys4};
                app.whchS=[2 4];
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==1 && includeS(4)==1
                Sys={Sys3,Sys4};
                app.whchS=[3 4];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==1 && includeS(4)==0
                Sys={Sys1,Sys2,Sys3};
                app.whchS=[1 2 3];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==0 && includeS(4)==1
                Sys={Sys1,Sys2,Sys4};
                app.whchS=[1 2 4];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==1 && includeS(4)==1
                Sys={Sys1,Sys3,Sys4};
                app.whchS=[1 3 4];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==1 && includeS(4)==1
                Sys={Sys2,Sys3,Sys4};
                app.whchS=[2 3 4];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==1 && includeS(4)==1
                Sys={Sys1,Sys2,Sys3,Sys4};
                app.whchS=[1 2 3 4];
            end
            
            for i=1:length(Sys)
                if isempty(fncmeth{i,1})
                    Opt.fnc{i}=funcval;
                    Opt.meth{i}=methodval;
                else
                    Opt.fnc{i}=fncmeth{i,1};
                    Opt.meth{i}=fncmeth{i,2};
                end
            end
            
            switch app.dataexists
                case 'y'
                    if get(app.DataOnlyButton,'Value')==1
                        set(app.PlotBothButton,'Value',1)
                        set(app.DataOnlyButton,'Value',0)
                    end
                    if Exp.('nPoints')~=length(app.data.B)
                        Exp.nPoints=length(app.data.B);
                    end
                case 'n'
                    set(app.PlotBothButton,'Value',0)
                    set(app.DataOnlyButton,'Value',0)
                    set(app.SimOnlyButton,'Value',1)
            end
            
            [app.data.simB,sim_spc]=app.custom_twomethods(Sys,Exp,Opt);
            app.simexists='y';
            app.tags.noisysim='n';
            app.tags.lastprocess='simulation';
            
            app.data.simspc=zeros(Exp.nPoints,4);
            k=1;
            for i=1:4
                switch onORoff{i}
                    case 'On'
                        app.data.simspc(:,i)=sim_spc(:,k);
                        k=k+1;
                end
            end
            
            showselect=get(app.ShowButtonGroup,'SelectedObject');
            selbutton=showselect.Text;
            showplots(app,selbutton,app.whchS,onORoff)
            
            app.data.orgsimspc=app.data.simspc;
            app.data.orgsimB=app.data.simB;
            set(app.gvaluesButton,'Value',0)
            
            tm=num2str(toc,'%.2f');
            set(app.SimTimeLabel,'Text',['sim time =',tm,' s'],'Visible','On')
            
        end

        % Value changed function: FunctionDropDown
        function FunctionSelection(app, event)
            switch app.FunctionDropDown.Value
                case 'garlic'
                    set(app.SimMethodDropDown,'Items',{'exact','perturb','perturb1',...
                        'perturb2','perturb3','perturb4','perturb5'},'Value','perturb2')
                case 'chili'
                    set(app.SimMethodDropDown,'Items',{'general','Freed','Fast'},'Value','general')
                case 'pepper'
                    set(app.SimMethodDropDown,'Items',{'matrix','hybrid','perturb','perturb1','perturb2'},'Value','matrix')
            end
            
        end

        % Button pushed function: SaveParametersButton
        function saveparams(app, event)
            pars.Syslst1=get(app.Sys1ListBox,'Items');
            pars.Syslst2=get(app.Sys2ListBox,'Items');
            pars.Syslst3=get(app.Sys3ListBox,'Items');
            pars.Syslst4=get(app.Sys4ListBox,'Items');
            
            pars.Explst=get(app.ExpListBox,'Items');
            pars.Optlst=get(app.OptListBox,'Items');
            pars.FitOptlst=get(app.FitOptListBox,'Items');
            pars.Vary1lst=get(app.Vary1ListBox,'Items');
            pars.Vary2lst=get(app.Vary2ListBox,'Items');
            pars.Vary3lst=get(app.Vary3ListBox,'Items');
            pars.Vary4lst=get(app.Vary4ListBox,'Items');
            
            pars.funcitems=get(app.FunctionDropDown,'Items');
            pars.funcval=get(app.FunctionDropDown,'Value');
            
            pars.methitems=get(app.SimMethodDropDown,'Items');
            pars.methval=get(app.SimMethodDropDown,'Value');
            
            app.cwEPRapp.Visible = 'Off';
            try
                uisave('pars')
            catch
                app.cwEPRapp.Visible = 'On';
                return
            end
            app.cwEPRapp.Visible = 'On';
            
        end

        % Button pushed function: LoadParametersButton
        function loadparams(app, event)
            app.cwEPRapp.Visible = 'Off';
            try
                [name,path]=uigetfile('*.mat','Multiselect','Off');
                in=load(strcat(path,filesep,name));
            catch
                app.cwEPRapp.Visible = 'On';
                return
            end
            app.cwEPRapp.Visible = 'On';
            
            set(app.Sys1ListBox,'Items',in.pars.Syslst1);
            app.inp.Sys1Str=in.pars.Syslst1;
            set(app.Sys2ListBox,'Items',in.pars.Syslst2);
            app.inp.Sys2Str=in.pars.Syslst2;
            set(app.Sys3ListBox,'Items',in.pars.Syslst3);
            app.inp.Sys3Str=in.pars.Syslst3;
            set(app.Sys4ListBox,'Items',in.pars.Syslst4);
            app.inp.Sys4Str=in.pars.Syslst4;
            
            set(app.ExpListBox,'Items',in.pars.Explst);
            app.inp.ExpStr=in.pars.Explst;
            set(app.OptListBox,'Items',in.pars.Optlst);
            app.inp.OptStr=in.pars.Optlst;
            set(app.FitOptListBox,'Items',in.pars.FitOptlst);
            app.inp.FitOptStr=in.pars.FitOptlst;
            set(app.Vary1ListBox,'Items',in.pars.Vary1lst);
            app.inp.Vary1Str=in.pars.Vary1lst;
            set(app.Vary2ListBox,'Items',in.pars.Vary2lst);
            app.inp.Vary2Str=in.pars.Vary2lst;
            set(app.Vary3ListBox,'Items',in.pars.Vary3lst);
            app.inp.Vary3Str=in.pars.Vary3lst;
            set(app.Vary4ListBox,'Items',in.pars.Vary4lst);
            app.inp.Vary4Str=in.pars.Vary4lst;
            
            set(app.FunctionDropDown,'Items',in.pars.funcitems)
            set(app.FunctionDropDown,'Value',in.pars.funcval)
            
            set(app.SimMethodDropDown,'Items',in.pars.methitems)
            set(app.SimMethodDropDown,'Value',in.pars.methval)
            
            
        end

        % Button pushed function: esfitButton
        function esfit(app, event)
            switch app.dataexists
                case 'n'
                    set(app.FittingSliceLabel,"Visible","On","Text","There is no data to fit.")
                    return
            end
            onORoffS{1}=get(app.Sys1Switch,'Value');
            onORoffS{2}=get(app.Sys2Switch,'Value');
            onORoffS{3}=get(app.Sys3Switch,'Value');
            onORoffS{4}=get(app.Sys4Switch,'Value');
            
            onORoffV{1}=get(app.Vary1Switch,'Value');
            onORoffV{2}=get(app.Vary2Switch,'Value');
            onORoffV{3}=get(app.Vary3Switch,'Value');
            onORoffV{4}=get(app.Vary4Switch,'Value');
            
            Vary={};
            Sys={};
            s=0;
            switch onORoffS{1}
                case 'On'
                    syslst1=get(app.Sys1ListBox,'Items');
                    s=s+1;
                    [Sys1,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst1);
                    if isempty(Sys1)
                        set(app.Sys1Switch,'Value','Off')
                        Sys1=[];
                        includeS(1)=0;
                    else
                        includeS(1)=1;
                    end
                    switch onORoffV{1}
                        case 'On'
                            varypars_1=get(app.Vary1ListBox,'Items');
                            if isempty(varypars_1)
                                set(app.Vary1Switch,'Value','Off')
                                Vary1=[];
                            elseif ~iscell(varypars_1)
                                t=1;
                                varypars1{1}=varypars_1;
                                Vary1=getvarys(app,varypars1,t);
                            else
                                t=numel(varypars_1);
                                Vary1=getvarys(app,varypars_1,t);
                            end
                        case 'Off'
                            Vary1=[];
                    end
                case 'Off'
                    Sys1=[];
                    includeS(1)=0;
            end
            
            switch onORoffS{2}
                case 'On'
                    syslst2=get(app.Sys2ListBox,'Items');
                    s=s+1;
                    [Sys2,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst2);
                    if isempty(Sys2)
                        set(app.Sys2Switch,'Value','Off')
                        Sys2=[];
                        includeS(2)=0;
                    else
                        includeS(2)=1;
                    end
                    switch onORoffV{2}
                        case 'On'
                            varypars_2=get(app.Vary2ListBox,'Items');
                            if isempty(varypars_2)
                                set(app.Vary2Switch,'Value','Off')
                                Vary2=[];
                            elseif ~iscell(varypars_2)
                                t=1;
                                varypars2{1}=varypars_2;
                                Vary2=getvarys(app,varypars2,t);
                            else
                                t=numel(varypars_2);
                                Vary2=getvarys(app,varypars_2,t);
                            end
                        case 'Off'
                            Vary2=[];
                    end
                case 'Off'
                    Sys2=[];
                    includeS(2)=0;
            end
            
            switch onORoffS{3}
                case 'On'
                    syslst3=get(app.Sys3ListBox,'Items');
                    s=s+1;
                    [Sys3,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst3);
                    if isempty(Sys3)
                        set(app.Sys3Switch,'Value','Off')
                        Sys3=[];
                        includeS(3)=0;
                    else
                        includeS(3)=1;
                    end
                    switch onORoffV{3}
                        case 'On'
                            varypars_3=get(app.Vary3ListBox,'Items');
                            if isempty(varypars_3)
                                set(app.Vary3Switch,'Value','Off')
                                Vary3=[];
                            elseif ~iscell(varypars_3)
                                t=1;
                                varypars3{1}=varypars_3;
                                Vary3=getvarys(app,varypars3,t);
                            else
                                t=numel(varypars_3);
                                Vary3=getvarys(app,varypars_3,t);
                            end
                        case 'Off'
                            Vary3=[];
                    end
                case 'Off'
                    Sys3=[];
                    includeS(3)=0;
            end
            
            switch onORoffS{4}
                case 'On'
                    syslst4=get(app.Sys4ListBox,'Items');
                    s=s+1;
                    [Sys4,fncmeth{s,1},fncmeth{s,2}]=getsyss(app,syslst4);
                    if isempty(Sys4)
                        set(app.Sys4Switch,'Value','Off')
                        Sys4=[];
                        includeS(4)=0;
                    else
                        includeS(4)=1;
                    end
                    switch onORoffV{4}
                        case 'On'
                            varypars_4=get(app.Vary4ListBox,'Items');
                            if isempty(varypars_4)
                                set(app.Vary4Switch,'Value','Off')
                                Vary4=[];
                            elseif ~iscell(varypars_4)
                                t=1;
                                varypars4{1}=varypars_4;
                                Vary4=getvarys(app,varypars4,t);
                            else
                                t=numel(varypars_4);
                                Vary4=getvarys(app,varypars_4,t);
                            end
                        case 'Off'
                            Vary4=[];
                    end
                case 'Off'
                    Sys4=[];
                    includeS(4)=0;
            end
            
            if sum(includeS)==0
                set(app.FittingSliceLabel,'Visible','On','Text','No spin systems')
                return
            else
                set(app.FittingSliceLabel,'Visible','Off')
            end
            
            explst=get(app.ExpListBox,'Items');
            Exp=getexps(app,explst);
            
            if ~isfield(Exp,'mwFreq')
                set(app.SimTimeLabel,'Text','give mwFreq')
                return
            end
            app.params.simFq=Exp.('mwFreq');
            
            optlst=get(app.OptListBox,'Items');
            Opt = getopts(app,optlst);
            
            fitoptlst=get(app.FitOptListBox,'Items');
            [FitOpt, fitrange] = getfitopts(app,fitoptlst);
            
            
            if includeS(1)==1 && includeS(2)==0 && includeS(3)==0 && includeS(4)==0
                Sys=Sys1;
                Vary=Vary1;
                app.whchV=1;
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==0 && includeS(4)==0
                Sys=Sys2;
                Vary=Vary2;
                app.whchV=2;
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==1 && includeS(4)==0
                Sys=Sys3;
                Vary=Vary3;
                app.whchV=3;
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==0 && includeS(4)==1
                Sys=Sys4;
                Vary=Vary4;
                app.whchV=4;
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==0 && includeS(4)==0
                Sys={Sys1,Sys2};
                Vary={Vary1,Vary2};
                app.whchV=[1 2];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==1 && includeS(4)==0
                Sys={Sys1,Sys3};
                Vary={Vary1,Vary3};
                app.whchV=[1 3];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==0 && includeS(4)==1
                Sys={Sys1,Sys4};
                Vary={Vary1,Vary4};
                app.whchV=[1 4];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==1 && includeS(4)==0
                Sys={Sys2,Sys3};
                Vary={Vary2,Vary3};
                app.whchV=[2 3];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==0 && includeS(4)==1
                Sys={Sys2,Sys4};
                Vary={Vary2,Vary4};
                app.whchV=[2 4];
            elseif includeS(1)==0 && includeS(2)==0 && includeS(3)==1 && includeS(4)==1
                Sys={Sys3,Sys4};
                Vary={Vary3,Vary4};
                app.whchV=[3 4];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==1 && includeS(4)==0
                Sys={Sys1,Sys2,Sys3};
                Vary={Vary1,Vary2,Vary3};
                app.whchV=[1 2 3];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==0 && includeS(4)==1
                Sys={Sys1,Sys2,Sys4};
                Vary={Vary1,Vary2,Vary4};
                app.whchV=[1 2 4];
            elseif includeS(1)==1 && includeS(2)==0 && includeS(3)==1 && includeS(4)==1
                Sys={Sys1,Sys3,Sys4};
                Vary={Vary1,Vary3,Vary4};
                app.whchV=[1 3 4];
            elseif includeS(1)==0 && includeS(2)==1 && includeS(3)==1 && includeS(4)==1
                Sys={Sys2,Sys3,Sys4};
                Vary={Vary2,Vary3,Vary4};
                app.whchV=[2 3 4];
            elseif includeS(1)==1 && includeS(2)==1 && includeS(3)==1 && includeS(4)==1
                Sys={Sys1,Sys2,Sys3,Sys4};
                Vary={Vary1,Vary2,Vary3,Vary4};
                app.whchV=[1 2 3 4];
            end
            
            ftoptMethodval=get(app.FitMethodDropDown,'Value');
            ftopttargetval=get(app.TargetDropDown,'Value');
            
            switch ftopttargetval
                case 'as is'
                    mth='fcn';
                case 'integral'
                    mth='int';
                case 'dbl integral'
                    mth='dint';
                case 'derivative'
                    mth='diff';
                case 'FFT'
                    mth='fft';
            end
            FitOpt.Method = [ftoptMethodval,' ',mth];
            
            FitOpt.Scaling=get(app.FitScalingDropDown,'Value');
            
            fitfuncval=get(app.FunctionDropDown,'Value');
            fitmethodval=get(app.SimMethodDropDown,'Value');
            
            for i=1:s
                if isempty(fncmeth{i,1})
                    Opt.fnc{i}=fitfuncval;
                    Opt.meth{i}=fitmethodval;
                else
                    Opt.fnc{i}=fncmeth{i,1};
                    Opt.meth{i}=fncmeth{i,2};
                end
            end
            
            xval=app.data.B;
            spec=app.data.spc;
            
            if ~isempty(fitrange)
                rangeL=find(round(app.data.B)==fitrange(1));
                Opt.rangeLH(1)=rangeL(1);
                rangeH=find(round(app.data.B)==fitrange(2));
                Opt.rangeLH(2)=rangeH(end);
                for i=1:length(spec(1,:))
                    spc(:,i)=zeros(length(spec(:,1)),1);
                    spc(Opt.rangeLH(1):Opt.rangeLH(2),i)=spec(Opt.rangeLH(1):Opt.rangeLH(2),i);
                end
            else
                spc=spec;
                Opt.rangeLH=[1 length(app.data.B)];
            end
            
            if length(app.data.spc(1,:))==1
                fitselection='Only Current';
            else
                fitselect=get(app.DDataButtonGroup,'SelectedObject');
                fitselection=fitselect.Text;
            end
            
            v=0;
            if iscell(Vary)
                for i=1:length(Vary)
                    if isempty(Vary{i})
                        v=v+1;
                    end
                end
            end
            if isempty(Vary) || v==s
                set(app.FittingSliceLabel,'Visible','On','Text','Nothing Varying')
                return
            else
                set(app.FittingSliceLabel,'Visible','Off')
            end
            
            if Exp.('nPoints')~=length(app.data.B)
                Exp.nPoints=length(app.data.B);
            end
            switch fitselection
                case 'Only Current'
                    switch get(app.FitGUICheckBox,'Value')
                        case 1
                            esfit(@app.custom_fittwomethods,spc(:,app.currslice),Sys,Vary,Exp,Opt,FitOpt)
                        case 0
                            set(app.FittingSliceLabel,'Visible','On','Text',['Fitting Current Slice: ',num2str(app.currslice),' (see command window)'])
                            drawnow
                            FitOpt.PrintLevel=1;
                            [pars,app.data.fitspc,~]=esfit(@app.custom_fittwomethods,spc(:,app.currslice),Sys,Vary,Exp,Opt,FitOpt);
                            disp(pars)
                            fit1.Sys=pars;
                            assignin('base','fit1',fit1)
                            set(app.FittingSliceLabel,'Text','Fitting finished')
                            plot(app.UIAxes,xval,spc(:,app.currslice),'linewidth',app.dlw,'color',app.dclr);
                            hold(app.UIAxes,'on')
                            plot(app.UIAxes,xval,app.data.fitspc,'linewidth',app.slw,'color',app.sclr{1});
                            xlim(app.UIAxes,[min(xval) max(xval)]);
                            ylim(app.UIAxes,[min(spc(:,app.currslice)) max(spc(:,app.currslice))].*1.05);
                            legend(app.UIAxes,'Data','Sim','location','northeast');...
                                xlabel(app.UIAxes,app.xlab);ylabel(app.UIAxes,'Amplitude (a.u.)');
                            hold(app.UIAxes,'off')
                    end
                case 'All Sequential'
                    set(app.StopFittingButton,'Visible','On')
                    set(app.FitGUICheckBox,'Value',0)
                    FitOpt.PrintLevel=0;
                    kk=length(spc(1,:));
                    xlabel(app.UIAxes,app.xlab);
                    ylabel(app.UIAxes,'Amplitude (a.u.)');
                    set(app.FittingSliceLabel,'Visible','On')
                    app.data.fitspc=zeros(length(spc(:,1)),kk);
                    for k=1:kk
                        if get(app.StopFittingButton,'Value')==1
                            set(app.StopFittingButton,'Visible','Off','Value',0)
                            set(app.FittingSliceLabel,'Text','Fitting stopped')
                            return
                        end
                        set(app.FittingSliceLabel,'Text',['Fitting Slice: ',num2str(k),' of ',num2str(kk)])
                        set(app.Spinner,'Value',k)
                        drawnow
                        [pars,app.data.fitspc(:,k),~]=esfit(@app.custom_fittwomethods,spc(:,k),Sys,Vary,Exp,Opt,FitOpt);
                        disp(['Slice ',num2str(k),':'])
                        disp(pars)
                        fit1=pars;
                        assignin('base',['fit',num2str(k)],fit1)
                        plot(app.UIAxes,xval,spc(:,k),'linewidth',app.dlw,'color',app.dclr);
                        hold(app.UIAxes,'on')
                        plot(app.UIAxes,xval,app.data.fitspc(:,k),'linewidth',app.slw,'color',app.sclr{1});...
                            xlim(app.UIAxes,[min(xval) max(xval)]);
                        ylim(app.UIAxes,[min(spc(:,k)) max(spc(:,k))].*1.05);
                        legend(app.UIAxes,'Data','Fit','location','northeast');
                        hold(app.UIAxes,'off')
                        drawnow
                    end
                    set(app.StopFittingButton,'Visible','Off','Value',0)
                    set(app.FittingSliceLabel,'Text','Fitting finished')
                case 'All Parallel'
                    for i=1:length(spc(1,:))
                        k=((i-1)*length(spc(:,1)))+1;
                        kk=((i)*length(spc(:,1)));
                        EXspc(1,k:kk)=spc(:,i);
                    end
                    Exp.nSpcs=length(app.data.spc(1,:));
                    Exp.nP=Exp.nPoints;
                    Exp.nPoints=Exp.nSpcs*Exp.nP;
                    switch get(app.FitGUICheckBox,'Value')
                        case 1
                            set(app.FittingSliceLabel,'Text',' ','Visible','Off')
                            esfit(@app.custom_fit2D,EXspc,Sys,Vary,Exp,Opt,FitOpt);
                        case 0
                            set(app.FittingSliceLabel,'Text','Fitting all in parallel (see command window)','Visible','On')
                            drawnow
                            FitOpt.PrintLevel=1;
                            [pars,fit_spc,~]=esfit(@app.custom_fit2D,EXspc,Sys,Vary,Exp,Opt,FitOpt);
                            disp(pars)
                            fit1=pars;
                            assignin('base','fit1',fit1)
                            set(app.FittingSliceLabel,'Text','FITTING finished')
                            l=((app.currslice-1)*length(spc(:,1)))+1;
                            ll=((app.currslice)*length(spc(:,1)));
                            app.data.fitspc=fit_spc(l:ll);
                            plot(app.UIAxes,xval,EXspc(l:ll),xval,app.data.fitspc,'linewidth',app.slw,'color',app.sclr{1});...
                                xlim(app.UIAxes,[min(xval) max(xval)]);
                            ylim(app.UIAxes,[min(spc(:)) max(spc(:))].*1.05);
                            legend(app.UIAxes,'Data','Fit','location','northeast');...
                                xlabel(app.UIAxes,app.xlab);
                            ylabel(app.UIAxes,'Amplitude (a.u.)');
                    end
            end
            app.fitexists='y';
            
        end

        % Button pushed function: getesfitfit1Button
        function getesfitresults(app, event)
            
            sysstrct=evalin('base','fit1.Sys');
            
            nmsys=numel(sysstrct);
            tt=1;
            for t=app.whchV
                if t==1
                    syslst=get(app.Sys1ListBox,'Items');
                elseif t==2
                    syslst=get(app.Sys2ListBox,'Items');
                elseif t==3
                    syslst=get(app.Sys3ListBox,'Items');
                elseif t==4
                    syslst=get(app.Sys4ListBox,'Items');
                end
                
                if nmsys==1
                    sys_strct=sysstrct;
                else
                    sys_strct=sysstrct{tt};
                end
                sys_params=fieldnames(sys_strct);
                
                [syslst_new] = getsysstructs(app,syslst,sys_strct,sys_params);
                
                if t==1
                    app.inp.Sys1=syslst_new;
                    set(app.Sys1Edit,'Value','')
                    set(app.Sys1ListBox,'Items',syslst_new)
                elseif t==2
                    app.inp.Sys2=syslst_new;
                    set(app.Sys2Edit,'Value','')
                    set(app.Sys2ListBox,'Items',syslst_new)
                elseif t==3
                    app.inp.Sys3=syslst_new;
                    set(app.Sys3Edit,'Value','')
                    set(app.Sys3ListBox,'Items',syslst_new)
                elseif t==4
                    app.inp.Sys4=syslst_new;
                    set(app.Sys4Edit,'Value','')
                    set(app.Sys4ListBox,'Items',syslst_new)
                end
                tt=tt+1;
            end
            
        end

        % Value changed function: StopFittingButton
        function StopFitting(app, event)
            switch app.StopFittingButton.Value
                case 1
                    set(app.FittingSliceLabel,'Text',['Stopping after slice: ',num2str(app.currslice)])
                    drawnow
            end
            
        end

        % Callback function: RescaleDropDown, ShowButtonGroup
        function showDatas(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            selectedButton = app.ShowButtonGroup.SelectedObject;
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
            
        end

        % Button pushed function: ExportDataButton
        function Export(app, event)
            currdir=cd;
            show=get(app.ShowButtonGroup,'SelectedObject');
            
            switch show.Text
                case 'Sim Only'
                    switch app.simexists
                        case 'y'
                            xval=app.data.simB;
                            yval=sum(app.data.simspc,2);
                            assignin('base','yval',yval)
                            mwFq=str2double(app.params.simFq);
                        case 'n'
                            xval=app.data.B;
                            yval=app.data.spc;
                            mwFq=str2double(app.params.MFQ);
                    end
                otherwise
                    switch app.dataexists
                        case 'n'
                            return
                    end
                    xval=app.data.B;
                    yval=app.data.spc;
                    mwFq=str2double(app.params.MFQ);
            end
            
            switch get(app.ExportDTADSCCheckBox,'Value')
                case 1
                    app.cwEPRapp.Visible = 'Off';
                    try
                        [fname,pth] = uiputfile({'*.*'});
                    catch
                        app.cwEPRapp.Visible = 'On';
                        return
                    end
                    app.cwEPRapp.Visible = 'On';
                    cd(pth)
                    if length(yval(1,:))>1
                        xval={xval,app.data.scnD};
                    end
                    if isempty(mwFq)
                        eprsave(fname,xval,yval,fname)
                    else
                        eprsave(fname,xval,yval,fname,mwFq);
                    end
                    cd(currdir)
            end
            switch get(app.ExportasciiCheckBox,'Value')
                case 1
                    svdata=zeros(length(xval),length(yval(1,:))+1);
                    svdata(:,1)=xval;
                    svdata(:,2:end)=yval;
                    currdir=cd;
                    app.cwEPRapp.Visible = 'Off';
                    try
                        [fnm,pth] = uiputfile({'*'});
                    catch
                        app.cwEPRapp.Visible = 'On';
                        return
                    end
                    app.cwEPRapp.Visible = 'On';
                    cd(pth)
                    fname=[fnm,'.dat'];
                    writematrix(svdata,fname);
                    cd(currdir)
            end
            switch get(app.ExporttoworkspaceCheckBox,'Value')
                case 1
                    switch app.dataexists
                        case 'y'
                            assignin('base','data_spc',app.data.spc)
                            if length(app.data.spc(1,:))>1
                                B={app.data.B,app.data.scnD};
                                assignin('base','data_B',B)
                            else
                                assignin('base','data_B',app.data.B)
                            end
                            
                    end
                    switch app.simexists
                        case 'y'
                            assignin('base','sim_spc',app.data.simspc)
                            assignin('base','sim_B',app.data.simB)
                    end
                    assignin('base','params',app.params)
            end
            
        end

        % Button pushed function: ExportStructsButton
        function export_simStructs(app, event)
            
            explst=get(app.ExpListBox,'Items');
            cwepr.Exp = getexps(app,explst);
            if isempty(cwepr.Exp)
                cwepr = rmfield(cwepr, "Exp");
            end
            
            optlst=get(app.OptListBox,'Items');
            cwepr.Opt = getopts(app,optlst);
            if isempty(cwepr.Opt)
                cwepr = rmfield(cwepr, "Opt");
            end
            
            sys_lsts = {get(app.Sys1ListBox,'Items'),...
                get(app.Sys2ListBox,'Items'),...
                get(app.Sys3ListBox,'Items'),...
                get(app.Sys4ListBox,'Items')};
            for k = 1:4
                SysValues = {};
                SysString = {};
                syslst=sys_lsts{k};
                l=1;
                for i=1:numel(syslst)
                    syslst_f=strsplit(syslst{i},'=');
                    syslst_f{2} = erase(syslst_f{2},"'");
                    if ~isempty(syslst_f{2})
                        SysString{l}=syslst_f{1};
                        SysValues{l}=str2num(syslst_f{2});
                        if isempty(SysValues{l})
                            SysValues{l}=syslst_f{2};
                        end
                        l=l+1;
                    else
                    end
                end
                if ~isempty(SysValues)
                    if k==1
                        cwepr.Sys1 = cell2struct(SysValues,SysString,2);
                    elseif k==2
                        cwepr.Sys2 = cell2struct(SysValues,SysString,2);
                    elseif k==3
                        cwepr.Sys3 = cell2struct(SysValues,SysString,2);
                    elseif k==4
                        cwepr.Sys4 = cell2struct(SysValues,SysString,2);
                    end
                end
            end
            
            vary_lsts={get(app.Vary1ListBox,'Items'),...
                get(app.Vary2ListBox,'Items'),...
                get(app.Vary3ListBox,'Items'),...
                get(app.Vary4ListBox,'Items')};
            for k=1:4
                if ~isempty(vary_lsts{k})
                    vary_lst = vary_lsts{k};
                    if ~iscell(vary_lst)
                        t=1;
                        vary_lst{1}=vary_lst;
                        vry=getvarys(app,vary_lst,t);
                    else
                        t=numel(vary_lst);
                        vry=getvarys(app,vary_lst,t);
                    end
                    if k==1
                        cwepr.Vary1 = vry;
                    elseif k==2
                        cwepr.Vary2 = vry;
                    elseif k==3
                        cwepr.Vary3 = vry;
                    elseif k==4
                        cwepr.Vary4 = vry;
                    end
                end
            end
            
            fitoptlst=get(app.FitOptListBox,'Items');
            [cwepr.FitOpt, ~] = getfitopts(app,fitoptlst);
            if isempty(cwepr.FitOpt)
                cwepr = rmfield(cwepr, "FitOpt");
            end
            
            if ~isempty(fieldnames(cwepr))
                assignin("base","cwepr",cwepr)
            end
            
            
        end

        % Value changed function: ExportasciiCheckBox
        function asciichk(app, event)
            switch app.ExportasciiCheckBox.Value
                case 1
                    set(app.ExportDTADSCCheckBox,'Value',0)
            end
            
        end

        % Value changed function: ExportDTADSCCheckBox
        function dtadscchk(app, event)
            switch app.ExportDTADSCCheckBox.Value
                case 1
                    set(app.ExportasciiCheckBox,'Value',0)
            end
            
        end

        % Value changed function: RealCheckBox
        function RealCheck(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            
            switch app.RealCheckBox.Value
                case 1
                    set(app.ImaginaryCheckBox,'Value',0)
                    app.data.spc = app.data.orgrealspc;
                    app.data.orgspc = app.data.orgrealspc;
                case 0
                    set(app.ImaginaryCheckBox,'Value',1)
                    app.data.spc = app.data.orgimagspc;
                    app.data.orgspc = app.data.orgimagspc;
            end
            
            selectedButton = app.ShowButtonGroup.SelectedObject;
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
            
        end

        % Value changed function: ImaginaryCheckBox
        function ImagCheck(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            
            switch app.ImaginaryCheckBox.Value
                case 1
                    set(app.RealCheckBox,'Value',0)
                    app.data.spc = app.data.orgimagspc;
                    app.data.orgspc = app.data.orgimagspc;
                case 0
                    set(app.RealCheckBox,'Value',1)
                    app.data.spc = app.data.orgrealspc;
                    app.data.orgspc = app.data.orgrealspc;
            end
            
            selectedButton = app.ShowButtonGroup.SelectedObject;
            selbutton=selectedButton.Text;
            onORoff{1}=get(app.Sys1Switch,'Value');
            onORoff{2}=get(app.Sys2Switch,'Value');
            onORoff{3}=get(app.Sys3Switch,'Value');
            onORoff{4}=get(app.Sys4Switch,'Value');
            togglesysplts(app,onORoff,selbutton)
            
        end

        % Button pushed function: FindconcentrationButton
        function spin_count(app, event)
            switch app.dataexists
                case 'n'
                    return
            end
            
            a = app.aEditField.Value;
            b = app.bEditField.Value;
            c = app.cEditField.Value;
            
            if a==0 && b==0 && c==0
                return
            end
            
            app.Integrate(app)
            
            fac1 = app.DIfactorEditField.Value;
            fac1 = strrep(fac1,',','*');
            if app.params.QVal
                if contains(fac1,'Q')==1
                    fac1 = strrep(fac1,'Q',app.params.QVal);
                end
            end
            if app.params.POW
                if contains(fac1,'Power')==1
                    fac1 = strrep(fac1,'Power',...
                        num2str(sqrt(str2double(app.params.POW))));
                end
            end
            if app.params.modA
                if contains(fac1,'ModAmp')==1
                    fac1 = strrep(fac1,'ModAmp',app.params.modA);
                end
            end
            if app.params.Temp
                if contains(fac1,'Temp')==1
                    fac1 = strrep(fac1,'Temp',app.params.Temp);
                end
            end
            
            fac = str2num(fac1);
            if isempty(fac) || fac==0
                set(app.SpinCEditField,'Value','Invalid Factor')
                return
            end
            
            for k = 1:length(app.data.spc(1,:))
                if isempty(app.aEditField) || a == 0
                    spins(:,k) = ((app.data.dblintgrl(:,k)/fac)-c) / b;
                else
                    d = c - (app.data.dblintgrl(:,k)/fac);
                    spins1(:,k) = ((-1*b)+sqrt((b^2)-(4*a*d))) / (2*a);
                    spins2(:,k) = ((-1*b)-sqrt((b^2)-(4*a*d))) / (2*a);
                end
            end
            
            if isempty(app.aEditField) || a == 0
                set(app.SpinCEditField,'Value',num2str(spins(:,app.currslice)))
                app.data.spinC = spins;
            else
                set(app.SpinCEditField,'Value',...
                    [num2str(spins1(:,app.currslice)),', ',num2str(spins2(:,app.currslice))])
                app.data.spinC = spins1;
            end
            
            app.tags.spincounted = 'y';
            
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create cwEPRapp
            app.cwEPRapp = uifigure;
            app.cwEPRapp.AutoResizeChildren = 'off';
            app.cwEPRapp.Position = [100 100 1200 720];
            app.cwEPRapp.Name = 'UI Figure';
            app.cwEPRapp.Resize = 'off';

            % Create UIAxes
            app.UIAxes = uiaxes(app.cwEPRapp);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, 'g')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.PlotBoxAspectRatio = [2.01388888888889 1 1];
            app.UIAxes.Position = [12 262 774 417];

            % Create OpenDataButton
            app.OpenDataButton = uibutton(app.cwEPRapp, 'push');
            app.OpenDataButton.ButtonPushedFcn = createCallbackFcn(app, @Opendata, true);
            app.OpenDataButton.Position = [18.5 682 83 22];
            app.OpenDataButton.Text = 'Open Data...';

            % Create SysTabGroup
            app.SysTabGroup = uitabgroup(app.cwEPRapp);
            app.SysTabGroup.AutoResizeChildren = 'off';
            app.SysTabGroup.Position = [844 74 348 583];

            % Create Sys1Tab
            app.Sys1Tab = uitab(app.SysTabGroup);
            app.Sys1Tab.AutoResizeChildren = 'off';
            app.Sys1Tab.Title = 'Sys1';

            % Create Sys1Switch
            app.Sys1Switch = uiswitch(app.Sys1Tab, 'slider');
            app.Sys1Switch.ValueChangedFcn = createCallbackFcn(app, @sys1selswitch, true);
            app.Sys1Switch.Position = [268 526 45 20];
            app.Sys1Switch.Value = 'On';

            % Create Sys1Edit
            app.Sys1Edit = uieditfield(app.Sys1Tab, 'text');
            app.Sys1Edit.ValueChangedFcn = createCallbackFcn(app, @EditSys1, true);
            app.Sys1Edit.Position = [9 525 226 22];

            % Create Sys1ListBox
            app.Sys1ListBox = uilistbox(app.Sys1Tab);
            app.Sys1ListBox.Items = {};
            app.Sys1ListBox.ValueChangedFcn = createCallbackFcn(app, @Sys1List, true);
            app.Sys1ListBox.Position = [10 8 328 505];
            app.Sys1ListBox.Value = {};

            % Create Sys2Tab
            app.Sys2Tab = uitab(app.SysTabGroup);
            app.Sys2Tab.AutoResizeChildren = 'off';
            app.Sys2Tab.Title = 'Sys2';

            % Create Sys2Switch
            app.Sys2Switch = uiswitch(app.Sys2Tab, 'slider');
            app.Sys2Switch.ValueChangedFcn = createCallbackFcn(app, @sys2selswitch, true);
            app.Sys2Switch.Position = [268 526 45 20];
            app.Sys2Switch.Value = 'On';

            % Create Sys2Edit
            app.Sys2Edit = uieditfield(app.Sys2Tab, 'text');
            app.Sys2Edit.ValueChangedFcn = createCallbackFcn(app, @EditSys2, true);
            app.Sys2Edit.Position = [9 525 226 22];

            % Create Sys2ListBox
            app.Sys2ListBox = uilistbox(app.Sys2Tab);
            app.Sys2ListBox.Items = {};
            app.Sys2ListBox.ValueChangedFcn = createCallbackFcn(app, @Sys2List, true);
            app.Sys2ListBox.Position = [10 8 328 505];
            app.Sys2ListBox.Value = {};

            % Create Sys3Tab
            app.Sys3Tab = uitab(app.SysTabGroup);
            app.Sys3Tab.AutoResizeChildren = 'off';
            app.Sys3Tab.Title = 'Sys3';

            % Create Sys3Switch
            app.Sys3Switch = uiswitch(app.Sys3Tab, 'slider');
            app.Sys3Switch.ValueChangedFcn = createCallbackFcn(app, @sys3selswitch, true);
            app.Sys3Switch.Position = [268 526 45 20];
            app.Sys3Switch.Value = 'On';

            % Create Sys3Edit
            app.Sys3Edit = uieditfield(app.Sys3Tab, 'text');
            app.Sys3Edit.ValueChangedFcn = createCallbackFcn(app, @EditSys3, true);
            app.Sys3Edit.Position = [9 525 226 22];

            % Create Sys3ListBox
            app.Sys3ListBox = uilistbox(app.Sys3Tab);
            app.Sys3ListBox.Items = {};
            app.Sys3ListBox.ValueChangedFcn = createCallbackFcn(app, @Sys3List, true);
            app.Sys3ListBox.Position = [10 8 328 505];
            app.Sys3ListBox.Value = {};

            % Create Sys4Tab
            app.Sys4Tab = uitab(app.SysTabGroup);
            app.Sys4Tab.AutoResizeChildren = 'off';
            app.Sys4Tab.Title = 'Sys4';

            % Create Sys4Switch
            app.Sys4Switch = uiswitch(app.Sys4Tab, 'slider');
            app.Sys4Switch.ValueChangedFcn = createCallbackFcn(app, @sys4selswitch, true);
            app.Sys4Switch.Position = [268 526 45 20];
            app.Sys4Switch.Value = 'On';

            % Create Sys4Edit
            app.Sys4Edit = uieditfield(app.Sys4Tab, 'text');
            app.Sys4Edit.ValueChangedFcn = createCallbackFcn(app, @EditSys4, true);
            app.Sys4Edit.Position = [9 525 226 22];

            % Create Sys4ListBox
            app.Sys4ListBox = uilistbox(app.Sys4Tab);
            app.Sys4ListBox.Items = {};
            app.Sys4ListBox.ValueChangedFcn = createCallbackFcn(app, @Sys4List, true);
            app.Sys4ListBox.Position = [10 8 328 505];
            app.Sys4ListBox.Value = {};

            % Create ExpTab
            app.ExpTab = uitab(app.SysTabGroup);
            app.ExpTab.AutoResizeChildren = 'off';
            app.ExpTab.Title = 'Exp';

            % Create ExpEdit
            app.ExpEdit = uieditfield(app.ExpTab, 'text');
            app.ExpEdit.ValueChangedFcn = createCallbackFcn(app, @EditExp, true);
            app.ExpEdit.Position = [9 525 226 22];

            % Create FromDataButton
            app.FromDataButton = uibutton(app.ExpTab, 'push');
            app.FromDataButton.ButtonPushedFcn = createCallbackFcn(app, @getexp, true);
            app.FromDataButton.Position = [241 525 100 22];
            app.FromDataButton.Text = 'From Data';

            % Create ExpListBox
            app.ExpListBox = uilistbox(app.ExpTab);
            app.ExpListBox.Items = {};
            app.ExpListBox.ValueChangedFcn = createCallbackFcn(app, @ExpList, true);
            app.ExpListBox.Position = [10 209 328 304];
            app.ExpListBox.Value = {};

            % Create FrequencySliderLabel
            app.FrequencySliderLabel = uilabel(app.ExpTab);
            app.FrequencySliderLabel.HorizontalAlignment = 'right';
            app.FrequencySliderLabel.Position = [7 177 62 22];
            app.FrequencySliderLabel.Text = 'Frequency';

            % Create FrequencySlider
            app.FrequencySlider = uislider(app.ExpTab);
            app.FrequencySlider.Limits = [9.1 9.9];
            app.FrequencySlider.MajorTicks = [];
            app.FrequencySlider.MajorTickLabels = {};
            app.FrequencySlider.ValueChangedFcn = createCallbackFcn(app, @freqChanged, true);
            app.FrequencySlider.ValueChangingFcn = createCallbackFcn(app, @freqChanging, true);
            app.FrequencySlider.MinorTicks = [];
            app.FrequencySlider.Position = [170 186 162 3];
            app.FrequencySlider.Value = 9.1;

            % Create StaticExpListBox
            app.StaticExpListBox = uilistbox(app.ExpTab);
            app.StaticExpListBox.Items = {};
            app.StaticExpListBox.Position = [13 11 325 152];
            app.StaticExpListBox.Value = {};

            % Create FreqSpinner
            app.FreqSpinner = uispinner(app.ExpTab);
            app.FreqSpinner.Step = 0.0001;
            app.FreqSpinner.ValueChangedFcn = createCallbackFcn(app, @freqChanged, true);
            app.FreqSpinner.Position = [76 177 84 22];

            % Create OptTab
            app.OptTab = uitab(app.SysTabGroup);
            app.OptTab.AutoResizeChildren = 'off';
            app.OptTab.Title = 'Opt';

            % Create OptEdit
            app.OptEdit = uieditfield(app.OptTab, 'text');
            app.OptEdit.ValueChangedFcn = createCallbackFcn(app, @EditOpt, true);
            app.OptEdit.Position = [9 525 226 22];

            % Create OptListBox
            app.OptListBox = uilistbox(app.OptTab);
            app.OptListBox.Items = {};
            app.OptListBox.ValueChangedFcn = createCallbackFcn(app, @OptList, true);
            app.OptListBox.Position = [10 8 328 505];
            app.OptListBox.Value = {};

            % Create ProcessingTabGroup
            app.ProcessingTabGroup = uitabgroup(app.cwEPRapp);
            app.ProcessingTabGroup.AutoResizeChildren = 'off';
            app.ProcessingTabGroup.Position = [11 8 485 221];

            % Create BaselineTab
            app.BaselineTab = uitab(app.ProcessingTabGroup);
            app.BaselineTab.AutoResizeChildren = 'off';
            app.BaselineTab.Title = 'Baseline';

            % Create LoadButton
            app.LoadButton = uibutton(app.BaselineTab, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBackground, true);
            app.LoadButton.Position = [117 161 58 22];
            app.LoadButton.Text = 'Load';

            % Create BackgroundFilenameLabel
            app.BackgroundFilenameLabel = uilabel(app.BaselineTab);
            app.BackgroundFilenameLabel.Position = [179 161 291 22];
            app.BackgroundFilenameLabel.Text = 'BackgroundFilename';

            % Create SpectrumUsageSliderLabel
            app.SpectrumUsageSliderLabel = uilabel(app.BaselineTab);
            app.SpectrumUsageSliderLabel.HorizontalAlignment = 'right';
            app.SpectrumUsageSliderLabel.Position = [10 88 96 22];
            app.SpectrumUsageSliderLabel.Text = 'Spectrum Usage';

            % Create SpectrumUsageSlider
            app.SpectrumUsageSlider = uislider(app.BaselineTab);
            app.SpectrumUsageSlider.Limits = [0 50];
            app.SpectrumUsageSlider.MajorTicks = [5 10 15 20 25 30 35 40 45 50];
            app.SpectrumUsageSlider.ValueChangedFcn = createCallbackFcn(app, @PolynomialFit, true);
            app.SpectrumUsageSlider.Position = [119 97 224 3];
            app.SpectrumUsageSlider.Value = 50;

            % Create PolynomialFitButton
            app.PolynomialFitButton = uibutton(app.BaselineTab, 'push');
            app.PolynomialFitButton.ButtonPushedFcn = createCallbackFcn(app, @PolynomialFit, true);
            app.PolynomialFitButton.Position = [13 118 131 22];
            app.PolynomialFitButton.Text = 'Polynomial Fit';

            % Create ApplyButton
            app.ApplyButton = uibutton(app.BaselineTab, 'push');
            app.ApplyButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyBaseline, true);
            app.ApplyButton.Position = [399 144 60 37];
            app.ApplyButton.Text = 'Apply';

            % Create BackgroundButton
            app.BackgroundButton = uibutton(app.BaselineTab, 'push');
            app.BackgroundButton.ButtonPushedFcn = createCallbackFcn(app, @Background, true);
            app.BackgroundButton.Position = [13 161 100 22];
            app.BackgroundButton.Text = 'Background';

            % Create RevertBaseButton
            app.RevertBaseButton = uibutton(app.BaselineTab, 'push');
            app.RevertBaseButton.ButtonPushedFcn = createCallbackFcn(app, @RevertBaseline, true);
            app.RevertBaseButton.Position = [399 88 60 22];
            app.RevertBaseButton.Text = 'Revert';

            % Create fromendsLabel
            app.fromendsLabel = uilabel(app.BaselineTab);
            app.fromendsLabel.Position = [21 70 81 22];
            app.fromendsLabel.Text = '(% from ends)';

            % Create msbackadjButton
            app.msbackadjButton = uibutton(app.BaselineTab, 'push');
            app.msbackadjButton.ButtonPushedFcn = createCallbackFcn(app, @msbackadjustment, true);
            app.msbackadjButton.Position = [13 21 131 22];
            app.msbackadjButton.Text = 'msbackadj';

            % Create RequiresBioinformaticsToolboxLabel
            app.RequiresBioinformaticsToolboxLabel = uilabel(app.BaselineTab);
            app.RequiresBioinformaticsToolboxLabel.FontAngle = 'italic';
            app.RequiresBioinformaticsToolboxLabel.Position = [150 21 184 22];
            app.RequiresBioinformaticsToolboxLabel.Text = '(Requires Bioinformatics Toolbox)';

            % Create NormalizeBaseCheckBox
            app.NormalizeBaseCheckBox = uicheckbox(app.BaselineTab);
            app.NormalizeBaseCheckBox.Text = 'Normalize';
            app.NormalizeBaseCheckBox.Position = [400 117 76 22];

            % Create OrderSpinnerLabel
            app.OrderSpinnerLabel = uilabel(app.BaselineTab);
            app.OrderSpinnerLabel.HorizontalAlignment = 'right';
            app.OrderSpinnerLabel.Position = [148 117 36 22];
            app.OrderSpinnerLabel.Text = 'Order';

            % Create PolyOrderSpinner
            app.PolyOrderSpinner = uispinner(app.BaselineTab);
            app.PolyOrderSpinner.Limits = [0 6];
            app.PolyOrderSpinner.ValueChangedFcn = createCallbackFcn(app, @PolynomialFit, true);
            app.PolyOrderSpinner.Position = [188 117 71 22];

            % Create SmoothTab
            app.SmoothTab = uitab(app.ProcessingTabGroup);
            app.SmoothTab.AutoResizeChildren = 'off';
            app.SmoothTab.Title = 'Smooth';

            % Create SmoothingOptionsButtonGroup
            app.SmoothingOptionsButtonGroup = uibuttongroup(app.SmoothTab);
            app.SmoothingOptionsButtonGroup.AutoResizeChildren = 'off';
            app.SmoothingOptionsButtonGroup.Title = 'Smoothing Options';
            app.SmoothingOptionsButtonGroup.Position = [13 9 384 174];

            % Create MovingAverageButton
            app.MovingAverageButton = uiradiobutton(app.SmoothingOptionsButtonGroup);
            app.MovingAverageButton.Text = 'Moving Average';
            app.MovingAverageButton.Position = [11 128 109 22];
            app.MovingAverageButton.Value = true;

            % Create BinomialButton
            app.BinomialButton = uiradiobutton(app.SmoothingOptionsButtonGroup);
            app.BinomialButton.Text = 'Binomial';
            app.BinomialButton.Position = [11 105 68 22];

            % Create FlatButton
            app.FlatButton = uiradiobutton(app.SmoothingOptionsButtonGroup);
            app.FlatButton.Text = 'Flat';
            app.FlatButton.Position = [11 83 42 22];

            % Create SavitskyGolayButton
            app.SavitskyGolayButton = uiradiobutton(app.SmoothingOptionsButtonGroup);
            app.SavitskyGolayButton.Text = 'Savitsky-Golay';
            app.SavitskyGolayButton.Position = [11 61 103 22];

            % Create wdenoiseButton
            app.wdenoiseButton = uiradiobutton(app.SmoothingOptionsButtonGroup);
            app.wdenoiseButton.Text = 'wdenoise';
            app.wdenoiseButton.Position = [11 29 73 22];

            % Create pdifEditFieldLabel
            app.pdifEditFieldLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.pdifEditFieldLabel.HorizontalAlignment = 'right';
            app.pdifEditFieldLabel.FontSize = 10;
            app.pdifEditFieldLabel.Position = [124 61 28 22];
            app.pdifEditFieldLabel.Text = 'p, dif';

            % Create pdifEditField
            app.pdifEditField = uieditfield(app.SmoothingOptionsButtonGroup, 'text');
            app.pdifEditField.Position = [158 61 110 22];

            % Create LevelLabel
            app.LevelLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.LevelLabel.HorizontalAlignment = 'right';
            app.LevelLabel.Position = [34 6 37 22];
            app.LevelLabel.Text = 'Level ';

            % Create wdOptionsEditField
            app.wdOptionsEditField = uieditfield(app.SmoothingOptionsButtonGroup, 'text');
            app.wdOptionsEditField.Position = [209 6 59 22];

            % Create RequiresWavletToolboxLabel
            app.RequiresWavletToolboxLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.RequiresWavletToolboxLabel.FontAngle = 'italic';
            app.RequiresWavletToolboxLabel.Position = [87 29 143 22];
            app.RequiresWavletToolboxLabel.Text = '(Requires Wavlet Toolbox)';

            % Create orfloorlog2nPointsLabel
            app.orfloorlog2nPointsLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.orfloorlog2nPointsLabel.FontSize = 11;
            app.orfloorlog2nPointsLabel.FontAngle = 'italic';
            app.orfloorlog2nPointsLabel.Position = [73 7 132 22];
            app.orfloorlog2nPointsLabel.Text = '(< or = floor(log2(nPoints))';

            % Create SmoothingBusyLabel
            app.SmoothingBusyLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.SmoothingBusyLabel.HorizontalAlignment = 'center';
            app.SmoothingBusyLabel.Visible = 'off';
            app.SmoothingBusyLabel.Position = [305 74 60 22];
            app.SmoothingBusyLabel.Text = 'Busy...';

            % Create StepsSpinnerLabel
            app.StepsSpinnerLabel = uilabel(app.SmoothingOptionsButtonGroup);
            app.StepsSpinnerLabel.HorizontalAlignment = 'right';
            app.StepsSpinnerLabel.Position = [250 99 37 22];
            app.StepsSpinnerLabel.Text = 'Steps';

            % Create StepsSpinner
            app.StepsSpinner = uispinner(app.SmoothingOptionsButtonGroup);
            app.StepsSpinner.Limits = [0 Inf];
            app.StepsSpinner.ValueChangedFcn = createCallbackFcn(app, @ShowSmoothing, true);
            app.StepsSpinner.Position = [291 99 76 22];
            app.StepsSpinner.Value = 1;

            % Create ShowButton
            app.ShowButton = uibutton(app.SmoothingOptionsButtonGroup, 'push');
            app.ShowButton.ButtonPushedFcn = createCallbackFcn(app, @ShowSmoothing, true);
            app.ShowButton.Position = [291 125 76 22];
            app.ShowButton.Text = 'Show';

            % Create SmoothingApplyButton
            app.SmoothingApplyButton = uibutton(app.SmoothTab, 'push');
            app.SmoothingApplyButton.ButtonPushedFcn = createCallbackFcn(app, @ApplySmoothing, true);
            app.SmoothingApplyButton.Position = [413 144 60 37];
            app.SmoothingApplyButton.Text = 'Apply';

            % Create SmoothingRevertButton
            app.SmoothingRevertButton = uibutton(app.SmoothTab, 'push');
            app.SmoothingRevertButton.ButtonPushedFcn = createCallbackFcn(app, @RevertSmoothing, true);
            app.SmoothingRevertButton.Position = [413 114 60 22];
            app.SmoothingRevertButton.Text = 'Revert';

            % Create IntegrateTab
            app.IntegrateTab = uitab(app.ProcessingTabGroup);
            app.IntegrateTab.AutoResizeChildren = 'off';
            app.IntegrateTab.Title = 'Integrate';

            % Create IntegrateButton
            app.IntegrateButton = uibutton(app.IntegrateTab, 'push');
            app.IntegrateButton.ButtonPushedFcn = createCallbackFcn(app, @Integrate, true);
            app.IntegrateButton.Position = [351 159 62 22];
            app.IntegrateButton.Text = 'Integrate';

            % Create ShowIntButtonGroup
            app.ShowIntButtonGroup = uibuttongroup(app.IntegrateTab);
            app.ShowIntButtonGroup.AutoResizeChildren = 'off';
            app.ShowIntButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ShowIntegrals, true);
            app.ShowIntButtonGroup.Title = 'Show';
            app.ShowIntButtonGroup.Position = [347 8 123 111];

            % Create SingleButton
            app.SingleButton = uiradiobutton(app.ShowIntButtonGroup);
            app.SingleButton.Text = 'Single';
            app.SingleButton.Position = [11 65 55 22];

            % Create DoubleButton
            app.DoubleButton = uiradiobutton(app.ShowIntButtonGroup);
            app.DoubleButton.Text = 'Double';
            app.DoubleButton.Position = [11 44 60 22];

            % Create BothButton
            app.BothButton = uiradiobutton(app.ShowIntButtonGroup);
            app.BothButton.Text = 'Both';
            app.BothButton.Position = [11 22 48 22];
            app.BothButton.Value = true;

            % Create NeitherButton
            app.NeitherButton = uiradiobutton(app.ShowIntButtonGroup);
            app.NeitherButton.Text = 'Neither';
            app.NeitherButton.Position = [11 0 61 22];

            % Create BaselineCorrButtonGroup
            app.BaselineCorrButtonGroup = uibuttongroup(app.IntegrateTab);
            app.BaselineCorrButtonGroup.AutoResizeChildren = 'off';
            app.BaselineCorrButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @Integrate, true);
            app.BaselineCorrButtonGroup.Title = 'Baseline Correction for the first integral';
            app.BaselineCorrButtonGroup.Position = [10 112 274 75];

            % Create NoneButton
            app.NoneButton = uiradiobutton(app.BaselineCorrButtonGroup);
            app.NoneButton.Text = 'None';
            app.NoneButton.Position = [11 29 51 22];
            app.NoneButton.Value = true;

            % Create PolynomialorderButton
            app.PolynomialorderButton = uiradiobutton(app.BaselineCorrButtonGroup);
            app.PolynomialorderButton.Text = 'Polynomial, order:';
            app.PolynomialorderButton.Position = [11 6 119 22];

            % Create IntOrderSpinner
            app.IntOrderSpinner = uispinner(app.BaselineCorrButtonGroup);
            app.IntOrderSpinner.Limits = [0 6];
            app.IntOrderSpinner.ValueChangedFcn = createCallbackFcn(app, @integrate_spinner, true);
            app.IntOrderSpinner.Position = [134 5 56 22];

            % Create IntegralTextAreaLabel
            app.IntegralTextAreaLabel = uilabel(app.IntegrateTab);
            app.IntegralTextAreaLabel.HorizontalAlignment = 'right';
            app.IntegralTextAreaLabel.Position = [297 123 45 22];
            app.IntegralTextAreaLabel.Text = 'Integral';

            % Create DIEditField
            app.DIEditField = uitextarea(app.IntegrateTab);
            app.DIEditField.Position = [349 122 118 25];

            % Create DoubleLabel
            app.DoubleLabel = uilabel(app.IntegrateTab);
            app.DoubleLabel.Position = [295 137 44 22];
            app.DoubleLabel.Text = 'Double';

            % Create ndMomentButton
            app.ndMomentButton = uibutton(app.IntegrateTab, 'push');
            app.ndMomentButton.ButtonPushedFcn = createCallbackFcn(app, @secondMoment, true);
            app.ndMomentButton.Position = [16 9 100 22];
            app.ndMomentButton.Text = '2nd Moment';

            % Create TrapzCheckBox
            app.TrapzCheckBox = uicheckbox(app.IntegrateTab);
            app.TrapzCheckBox.ValueChangedFcn = createCallbackFcn(app, @trapz_select, true);
            app.TrapzCheckBox.Text = 'Trapz';
            app.TrapzCheckBox.Position = [422 168 51 22];
            app.TrapzCheckBox.Value = true;

            % Create SumCheckBox
            app.SumCheckBox = uicheckbox(app.IntegrateTab);
            app.SumCheckBox.ValueChangedFcn = createCallbackFcn(app, @sum_select, true);
            app.SumCheckBox.Text = 'Sum';
            app.SumCheckBox.Position = [422 150 47 22];

            % Create SecondMomentEdit
            app.SecondMomentEdit = uieditfield(app.IntegrateTab, 'text');
            app.SecondMomentEdit.Position = [125 8 153 22];

            % Create IntSpectrumUsageSlider
            app.IntSpectrumUsageSlider = uislider(app.IntegrateTab);
            app.IntSpectrumUsageSlider.Limits = [0 50];
            app.IntSpectrumUsageSlider.MajorTicks = [5 10 15 20 25 30 35 40 45 50];
            app.IntSpectrumUsageSlider.ValueChangedFcn = createCallbackFcn(app, @Integrate, true);
            app.IntSpectrumUsageSlider.Position = [15 73 262 3];
            app.IntSpectrumUsageSlider.Value = 10;

            % Create SpectrumUsagefromendsLabel
            app.SpectrumUsagefromendsLabel = uilabel(app.IntegrateTab);
            app.SpectrumUsagefromendsLabel.Position = [58 83 175 22];
            app.SpectrumUsagefromendsLabel.Text = 'Spectrum Usage (% from ends)';

            % Create AutoButton
            app.AutoButton = uibutton(app.IntegrateTab, 'push');
            app.AutoButton.ButtonPushedFcn = createCallbackFcn(app, @auto_percent, true);
            app.AutoButton.Position = [237 83 38 21];
            app.AutoButton.Text = 'Auto';

            % Create AutoBusyLabel
            app.AutoBusyLabel = uilabel(app.IntegrateTab);
            app.AutoBusyLabel.Visible = 'off';
            app.AutoBusyLabel.Position = [284 81 42 22];
            app.AutoBusyLabel.Text = 'Busy...';

            % Create QuantifyTab
            app.QuantifyTab = uitab(app.ProcessingTabGroup);
            app.QuantifyTab.AutoResizeChildren = 'off';
            app.QuantifyTab.Title = 'Quantify';

            % Create Step1Label
            app.Step1Label = uilabel(app.QuantifyTab);
            app.Step1Label.Position = [10 168 310 22];
            app.Step1Label.Text = 'Step 1: Set up the integration details in the Integrate tab';

            % Create Step2Label
            app.Step2Label = uilabel(app.QuantifyTab);
            app.Step2Label.Position = [10 148 404 22];
            app.Step2Label.Text = 'Step 2: Enter the coefficients of a calibration curve fit where x are known ';

            % Create DIEditFieldLabel
            app.DIEditFieldLabel = uilabel(app.QuantifyTab);
            app.DIEditFieldLabel.HorizontalAlignment = 'right';
            app.DIEditFieldLabel.Position = [14 107 31 22];
            app.DIEditFieldLabel.Text = 'DI = ';

            % Create aEditField
            app.aEditField = uieditfield(app.QuantifyTab, 'numeric');
            app.aEditField.Position = [50 107 100 22];

            % Create x2Label
            app.x2Label = uilabel(app.QuantifyTab);
            app.x2Label.Position = [153 107 26 22];
            app.x2Label.Text = 'x^2';

            % Create Label
            app.Label = uilabel(app.QuantifyTab);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [166 107 25 22];
            app.Label.Text = '+';

            % Create bEditField
            app.bEditField = uieditfield(app.QuantifyTab, 'numeric');
            app.bEditField.Position = [200 107 100 22];

            % Create xLabel
            app.xLabel = uilabel(app.QuantifyTab);
            app.xLabel.Position = [303 107 25 22];
            app.xLabel.Text = 'x';

            % Create EditField_2Label
            app.EditField_2Label = uilabel(app.QuantifyTab);
            app.EditField_2Label.HorizontalAlignment = 'right';
            app.EditField_2Label.Position = [304 107 25 22];
            app.EditField_2Label.Text = '+';

            % Create cEditField
            app.cEditField = uieditfield(app.QuantifyTab, 'numeric');
            app.cEditField.Position = [338 107 100 22];

            % Create Step2Label_2
            app.Step2Label_2 = uilabel(app.QuantifyTab);
            app.Step2Label_2.Position = [52 132 379 22];
            app.Step2Label_2.Text = 'concentrations and DI are the double integrals of the known spectra:';

            % Create Step4Label
            app.Step4Label = uilabel(app.QuantifyTab);
            app.Step4Label.Position = [10 6 349 22];
            app.Step4Label.Text = 'Step 4: Calculate the spin concentration in the same units as x:';

            % Create FindconcentrationButton
            app.FindconcentrationButton = uibutton(app.QuantifyTab, 'push');
            app.FindconcentrationButton.ButtonPushedFcn = createCallbackFcn(app, @spin_count, true);
            app.FindconcentrationButton.Position = [351 32 121 22];
            app.FindconcentrationButton.Text = 'Find concentration';

            % Create Step3Label
            app.Step3Label = uilabel(app.QuantifyTab);
            app.Step3Label.Position = [10 82 319 22];
            app.Step3Label.Text = 'Step 3: Enter the factor used to scale the double integral: ';

            % Create SpinCEditField
            app.SpinCEditField = uieditfield(app.QuantifyTab, 'text');
            app.SpinCEditField.Position = [355 6 114 22];

            % Create DIfactorEditField
            app.DIfactorEditField = uieditfield(app.QuantifyTab, 'text');
            app.DIfactorEditField.HorizontalAlignment = 'center';
            app.DIfactorEditField.Position = [50 60 298 22];

            % Create hintsLabel
            app.hintsLabel = uilabel(app.QuantifyTab);
            app.hintsLabel.FontSize = 11;
            app.hintsLabel.FontAngle = 'italic';
            app.hintsLabel.Position = [51 40 287 22];
            app.hintsLabel.Text = 'hints: this can be a formula, e.g. sqrt(2)*5000*exp(-1/295)';

            % Create hint2Label
            app.hint2Label = uilabel(app.QuantifyTab);
            app.hint2Label.FontSize = 11;
            app.hint2Label.FontAngle = 'italic';
            app.hint2Label.Position = [79 26 214 22];
            app.hint2Label.Text = 'the double integral is divided by this factor';

            % Create AdjustTab
            app.AdjustTab = uitab(app.ProcessingTabGroup);
            app.AdjustTab.AutoResizeChildren = 'off';
            app.AdjustTab.Title = 'Adjust';

            % Create AddNoiseButton
            app.AddNoiseButton = uibutton(app.AdjustTab, 'push');
            app.AddNoiseButton.ButtonPushedFcn = createCallbackFcn(app, @AddNoise, true);
            app.AddNoiseButton.Position = [16 111 100 22];
            app.AddNoiseButton.Text = 'Add Noise';

            % Create RevertAdjustmentsButton
            app.RevertAdjustmentsButton = uibutton(app.AdjustTab, 'push');
            app.RevertAdjustmentsButton.ButtonPushedFcn = createCallbackFcn(app, @RevertAdjustments, true);
            app.RevertAdjustmentsButton.Position = [413 114 60 22];
            app.RevertAdjustmentsButton.Text = 'Revert';

            % Create NoiseModelButtonGroup
            app.NoiseModelButtonGroup = uibuttongroup(app.AdjustTab);
            app.NoiseModelButtonGroup.AutoResizeChildren = 'off';
            app.NoiseModelButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @AddNoise, true);
            app.NoiseModelButtonGroup.Title = 'Noise Model';
            app.NoiseModelButtonGroup.Position = [139 50 100 84];

            % Create fButton
            app.fButton = uiradiobutton(app.NoiseModelButtonGroup);
            app.fButton.Text = '1/f';
            app.fButton.Position = [11 41 36 22];
            app.fButton.Value = true;

            % Create UniformButton
            app.UniformButton = uiradiobutton(app.NoiseModelButtonGroup);
            app.UniformButton.Text = 'Uniform';
            app.UniformButton.Position = [11 22 65 22];

            % Create GaussianButton
            app.GaussianButton = uiradiobutton(app.NoiseModelButtonGroup);
            app.GaussianButton.Text = 'Gaussian';
            app.GaussianButton.Position = [11 2 72 22];

            % Create ApplyAdjustmentsButton
            app.ApplyAdjustmentsButton = uibutton(app.AdjustTab, 'push');
            app.ApplyAdjustmentsButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyAdjustments, true);
            app.ApplyAdjustmentsButton.Position = [413 144 60 37];
            app.ApplyAdjustmentsButton.Text = 'Apply';

            % Create InterpolateButton
            app.InterpolateButton = uibutton(app.AdjustTab, 'push');
            app.InterpolateButton.ButtonPushedFcn = createCallbackFcn(app, @Interpolate, true);
            app.InterpolateButton.Position = [16 154 100 22];
            app.InterpolateButton.Text = 'Interpolate';

            % Create NoisySimButton
            app.NoisySimButton = uibutton(app.AdjustTab, 'push');
            app.NoisySimButton.ButtonPushedFcn = createCallbackFcn(app, @NoisySim, true);
            app.NoisySimButton.Position = [16 9 100 22];
            app.NoisySimButton.Text = 'Noisy Sim';

            % Create NoiseDTADSCCheckBox
            app.NoiseDTADSCCheckBox = uicheckbox(app.AdjustTab);
            app.NoiseDTADSCCheckBox.ValueChangedFcn = createCallbackFcn(app, @NoiseDTADSC, true);
            app.NoiseDTADSCCheckBox.Text = 'DTA/DSC';
            app.NoiseDTADSCCheckBox.Position = [127 9 73 22];

            % Create NoiseasciiCheckBox
            app.NoiseasciiCheckBox = uicheckbox(app.AdjustTab);
            app.NoiseasciiCheckBox.ValueChangedFcn = createCallbackFcn(app, @Noisetoascii, true);
            app.NoiseasciiCheckBox.Text = '-ascii';
            app.NoiseasciiCheckBox.Position = [205 9 51 22];

            % Create NoisetoworkspaceCheckBox
            app.NoisetoworkspaceCheckBox = uicheckbox(app.AdjustTab);
            app.NoisetoworkspaceCheckBox.Text = 'to workspace';
            app.NoisetoworkspaceCheckBox.Position = [263 9 95 22];

            % Create NoisetodataCheckBox
            app.NoisetodataCheckBox = uicheckbox(app.AdjustTab);
            app.NoisetodataCheckBox.ValueChangedFcn = createCallbackFcn(app, @Noisetodata, true);
            app.NoisetodataCheckBox.Text = 'to data';
            app.NoisetodataCheckBox.Position = [361 9 60 22];

            % Create SNRSpinnerLabel
            app.SNRSpinnerLabel = uilabel(app.AdjustTab);
            app.SNRSpinnerLabel.HorizontalAlignment = 'right';
            app.SNRSpinnerLabel.Position = [7 77 30 22];
            app.SNRSpinnerLabel.Text = 'SNR';

            % Create SNRSpinner
            app.SNRSpinner = uispinner(app.AdjustTab);
            app.SNRSpinner.Limits = [1 Inf];
            app.SNRSpinner.ValueChangedFcn = createCallbackFcn(app, @AddNoise, true);
            app.SNRSpinner.Position = [44 77 75 22];
            app.SNRSpinner.Value = 100;

            % Create PointsSpinnerLabel
            app.PointsSpinnerLabel = uilabel(app.AdjustTab);
            app.PointsSpinnerLabel.HorizontalAlignment = 'right';
            app.PointsSpinnerLabel.Position = [126 154 39 22];
            app.PointsSpinnerLabel.Text = 'Points';

            % Create PointsSpinner
            app.PointsSpinner = uispinner(app.AdjustTab);
            app.PointsSpinner.Limits = [1 Inf];
            app.PointsSpinner.ValueDisplayFormat = '%.0f';
            app.PointsSpinner.ValueChangedFcn = createCallbackFcn(app, @Interpolate, true);
            app.PointsSpinner.Position = [171 154 93 22];
            app.PointsSpinner.Value = 1;

            % Create ToolsTab
            app.ToolsTab = uitab(app.ProcessingTabGroup);
            app.ToolsTab.AutoResizeChildren = 'off';
            app.ToolsTab.Title = 'Tools';

            % Create ExportDataButton
            app.ExportDataButton = uibutton(app.ToolsTab, 'push');
            app.ExportDataButton.ButtonPushedFcn = createCallbackFcn(app, @Export, true);
            app.ExportDataButton.Position = [9 9 105 22];
            app.ExportDataButton.Text = 'Export Data';

            % Create ConstantsDropDownLabel
            app.ConstantsDropDownLabel = uilabel(app.ToolsTab);
            app.ConstantsDropDownLabel.HorizontalAlignment = 'right';
            app.ConstantsDropDownLabel.Position = [237 165 60 22];
            app.ConstantsDropDownLabel.Text = 'Constants';

            % Create ConstantsDropDown
            app.ConstantsDropDown = uidropdown(app.ToolsTab);
            app.ConstantsDropDown.Items = {'Constants...', 'amu', 'angstrom', 'avogadro', 'barn', 'bmagn', 'bohrrad', 'boltzm', 'clight', 'echarge', 'emass', 'eps0', 'evolt', 'faraday', 'gfree', 'hartree', 'hbar', 'molgas', 'mu0', 'nmagn', 'nmass', 'planck', 'pmass', 'rydberg', 'pi'};
            app.ConstantsDropDown.ValueChangedFcn = createCallbackFcn(app, @Constants, true);
            app.ConstantsDropDown.Position = [312 165 152 22];
            app.ConstantsDropDown.Value = 'Constants...';

            % Create ConversionToolButton
            app.ConversionToolButton = uibutton(app.ToolsTab, 'push');
            app.ConversionToolButton.ButtonPushedFcn = createCallbackFcn(app, @Convert, true);
            app.ConversionToolButton.Position = [122 129 107 22];
            app.ConversionToolButton.Text = 'Conversion Tool';

            % Create PeriodicTableButton
            app.PeriodicTableButton = uibutton(app.ToolsTab, 'push');
            app.PeriodicTableButton.ButtonPushedFcn = createCallbackFcn(app, @PeriodicTable, true);
            app.PeriodicTableButton.Position = [122 161 107 22];
            app.PeriodicTableButton.Text = 'Periodic Table';

            % Create ExportDTADSCCheckBox
            app.ExportDTADSCCheckBox = uicheckbox(app.ToolsTab);
            app.ExportDTADSCCheckBox.ValueChangedFcn = createCallbackFcn(app, @dtadscchk, true);
            app.ExportDTADSCCheckBox.Text = 'DTA/DSC';
            app.ExportDTADSCCheckBox.Position = [119 9 73 22];

            % Create ExportasciiCheckBox
            app.ExportasciiCheckBox = uicheckbox(app.ToolsTab);
            app.ExportasciiCheckBox.ValueChangedFcn = createCallbackFcn(app, @asciichk, true);
            app.ExportasciiCheckBox.Text = '-ascii';
            app.ExportasciiCheckBox.Position = [197 9 51 22];
            app.ExportasciiCheckBox.Value = true;

            % Create ExporttoworkspaceCheckBox
            app.ExporttoworkspaceCheckBox = uicheckbox(app.ToolsTab);
            app.ExporttoworkspaceCheckBox.Text = 'to workspace';
            app.ExporttoworkspaceCheckBox.Position = [253 9 95 22];

            % Create ShowButtonGroup
            app.ShowButtonGroup = uibuttongroup(app.ToolsTab);
            app.ShowButtonGroup.AutoResizeChildren = 'off';
            app.ShowButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @showDatas, true);
            app.ShowButtonGroup.Title = 'Show';
            app.ShowButtonGroup.Position = [13 92 100 91];

            % Create DataOnlyButton
            app.DataOnlyButton = uiradiobutton(app.ShowButtonGroup);
            app.DataOnlyButton.Text = 'Data Only';
            app.DataOnlyButton.Position = [11 45 75 22];

            % Create SimOnlyButton
            app.SimOnlyButton = uiradiobutton(app.ShowButtonGroup);
            app.SimOnlyButton.Text = 'Sim Only';
            app.SimOnlyButton.Position = [11 23 70 22];

            % Create PlotBothButton
            app.PlotBothButton = uiradiobutton(app.ShowButtonGroup);
            app.PlotBothButton.Text = 'Both';
            app.PlotBothButton.Position = [11 1 48 22];
            app.PlotBothButton.Value = true;

            % Create ResetDataButton
            app.ResetDataButton = uibutton(app.ToolsTab, 'push');
            app.ResetDataButton.ButtonPushedFcn = createCallbackFcn(app, @ResetProcessing, true);
            app.ResetDataButton.Position = [372 9 105 22];
            app.ResetDataButton.Text = 'Reset Data';

            % Create ConstantLabel
            app.ConstantLabel = uilabel(app.ToolsTab);
            app.ConstantLabel.Position = [312 137 152 22];
            app.ConstantLabel.Text = '';

            % Create LoadParametersButton
            app.LoadParametersButton = uibutton(app.ToolsTab, 'push');
            app.LoadParametersButton.ButtonPushedFcn = createCallbackFcn(app, @loadparams, true);
            app.LoadParametersButton.Position = [251 52 105 22];
            app.LoadParametersButton.Text = 'Load Parameters';

            % Create SaveParametersButton
            app.SaveParametersButton = uibutton(app.ToolsTab, 'push');
            app.SaveParametersButton.ButtonPushedFcn = createCallbackFcn(app, @saveparams, true);
            app.SaveParametersButton.Position = [130 52 105 22];
            app.SaveParametersButton.Text = 'Save Parameters';

            % Create getesfitfit1Button
            app.getesfitfit1Button = uibutton(app.ToolsTab, 'push');
            app.getesfitfit1Button.ButtonPushedFcn = createCallbackFcn(app, @getesfitresults, true);
            app.getesfitfit1Button.Position = [372 52 105 22];
            app.getesfitfit1Button.Text = 'get esfit ''fit1''';

            % Create ExportStructsButton
            app.ExportStructsButton = uibutton(app.ToolsTab, 'push');
            app.ExportStructsButton.ButtonPushedFcn = createCallbackFcn(app, @export_simStructs, true);
            app.ExportStructsButton.Position = [9 52 105 22];
            app.ExportStructsButton.Text = 'Export Structs';

            % Create DTab
            app.DTab = uitab(app.ProcessingTabGroup);
            app.DTab.AutoResizeChildren = 'off';
            app.DTab.Title = '2D';

            % Create PausesEditFieldLabel
            app.PausesEditFieldLabel = uilabel(app.DTab);
            app.PausesEditFieldLabel.HorizontalAlignment = 'right';
            app.PausesEditFieldLabel.Position = [84 155 54 22];
            app.PausesEditFieldLabel.Text = 'Pause (s)';

            % Create PausesEditField
            app.PausesEditField = uieditfield(app.DTab, 'numeric');
            app.PausesEditField.HorizontalAlignment = 'center';
            app.PausesEditField.Position = [142 155 53 22];
            app.PausesEditField.Value = 0.25;

            % Create plambdanPreAvgdeltadirEditFieldLabel
            app.plambdanPreAvgdeltadirEditFieldLabel = uilabel(app.DTab);
            app.plambdanPreAvgdeltadirEditFieldLabel.HorizontalAlignment = 'right';
            app.plambdanPreAvgdeltadirEditFieldLabel.FontSize = 10;
            app.plambdanPreAvgdeltadirEditFieldLabel.Position = [15 58 138 22];
            app.plambdanPreAvgdeltadirEditFieldLabel.Text = 'p, lambda, nPreAvg, delta, dir';

            % Create ewrlsoptions
            app.ewrlsoptions = uieditfield(app.DTab, 'text');
            app.ewrlsoptions.Position = [168 58 300 22];

            % Create ewrlsButton
            app.ewrlsButton = uibutton(app.DTab, 'push');
            app.ewrlsButton.ButtonPushedFcn = createCallbackFcn(app, @ewrls, true);
            app.ewrlsButton.Position = [15 86 70 22];
            app.ewrlsButton.Text = 'ewrls';

            % Create ewrlsBusyLabel
            app.ewrlsBusyLabel = uilabel(app.DTab);
            app.ewrlsBusyLabel.Visible = 'off';
            app.ewrlsBusyLabel.Position = [95 81 58 22];
            app.ewrlsBusyLabel.Text = 'Busy...';

            % Create ApplyewrlsButton
            app.ApplyewrlsButton = uibutton(app.DTab, 'push');
            app.ApplyewrlsButton.ButtonPushedFcn = createCallbackFcn(app, @Applyewrls, true);
            app.ApplyewrlsButton.Position = [15 30 69 22];
            app.ApplyewrlsButton.Text = 'Apply';

            % Create RevertewrlsButton
            app.RevertewrlsButton = uibutton(app.DTab, 'push');
            app.RevertewrlsButton.ButtonPushedFcn = createCallbackFcn(app, @Revertewrls, true);
            app.RevertewrlsButton.Position = [97 30 72 22];
            app.RevertewrlsButton.Text = 'Revert';

            % Create PlayButton
            app.PlayButton = uibutton(app.DTab, 'state');
            app.PlayButton.ValueChangedFcn = createCallbackFcn(app, @Play2D, true);
            app.PlayButton.Text = 'Play';
            app.PlayButton.Position = [14 155 66 22];

            % Create CompareewrlstooriginalLabel
            app.CompareewrlstooriginalLabel = uilabel(app.DTab);
            app.CompareewrlstooriginalLabel.HorizontalAlignment = 'right';
            app.CompareewrlstooriginalLabel.Position = [253 25 143 22];
            app.CompareewrlstooriginalLabel.Text = 'Compare ewrls to original';

            % Create ewrlsSpinner
            app.ewrlsSpinner = uispinner(app.DTab);
            app.ewrlsSpinner.Limits = [1 Inf];
            app.ewrlsSpinner.ValueChangedFcn = createCallbackFcn(app, @ewrlsCompare, true);
            app.ewrlsSpinner.Position = [404 25 64 22];
            app.ewrlsSpinner.Value = 1;

            % Create ewrlsedLabel
            app.ewrlsedLabel = uilabel(app.DTab);
            app.ewrlsedLabel.FontSize = 10;
            app.ewrlsedLabel.FontAngle = 'italic';
            app.ewrlsedLabel.Visible = 'off';
            app.ewrlsedLabel.Position = [18 5 228 18];
            app.ewrlsedLabel.Text = 'Current data was averaged with ewrls';

            % Create totalslicesLabel
            app.totalslicesLabel = uilabel(app.DTab);
            app.totalslicesLabel.FontSize = 14;
            app.totalslicesLabel.Visible = 'off';
            app.totalslicesLabel.Position = [18 129 177 22];
            app.totalslicesLabel.Text = '# total slices';

            % Create mainLabel
            app.mainLabel = uilabel(app.cwEPRapp);
            app.mainLabel.FontSize = 22;
            app.mainLabel.Position = [908 679 220 29];
            app.mainLabel.Text = 'cwEPR (easyspin.org)';

            % Create VaryTabGroup
            app.VaryTabGroup = uitabgroup(app.cwEPRapp);
            app.VaryTabGroup.AutoResizeChildren = 'off';
            app.VaryTabGroup.Position = [507 8 329 221];

            % Create Vary1Tab
            app.Vary1Tab = uitab(app.VaryTabGroup);
            app.Vary1Tab.AutoResizeChildren = 'off';
            app.Vary1Tab.Title = 'Vary1';

            % Create Vary1Switch
            app.Vary1Switch = uiswitch(app.Vary1Tab, 'slider');
            app.Vary1Switch.Position = [253 166 45 20];
            app.Vary1Switch.Value = 'On';

            % Create Vary1Edit
            app.Vary1Edit = uieditfield(app.Vary1Tab, 'text');
            app.Vary1Edit.ValueChangedFcn = createCallbackFcn(app, @EditVary1, true);
            app.Vary1Edit.Position = [7 165 217 22];

            % Create Vary1ListBox
            app.Vary1ListBox = uilistbox(app.Vary1Tab);
            app.Vary1ListBox.Items = {};
            app.Vary1ListBox.ValueChangedFcn = createCallbackFcn(app, @Vary1List, true);
            app.Vary1ListBox.Position = [7 8 271 150];
            app.Vary1ListBox.Value = {};

            % Create Vary1RemoveButton
            app.Vary1RemoveButton = uibutton(app.Vary1Tab, 'push');
            app.Vary1RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @Vary1Remove, true);
            app.Vary1RemoveButton.FontSize = 8;
            app.Vary1RemoveButton.Position = [281 128 44 18];
            app.Vary1RemoveButton.Text = 'Remove';

            % Create Vary2Tab
            app.Vary2Tab = uitab(app.VaryTabGroup);
            app.Vary2Tab.AutoResizeChildren = 'off';
            app.Vary2Tab.Title = 'Vary2';

            % Create Vary2Switch
            app.Vary2Switch = uiswitch(app.Vary2Tab, 'slider');
            app.Vary2Switch.Position = [253 166 45 20];
            app.Vary2Switch.Value = 'On';

            % Create Vary2Edit
            app.Vary2Edit = uieditfield(app.Vary2Tab, 'text');
            app.Vary2Edit.ValueChangedFcn = createCallbackFcn(app, @EditVary2, true);
            app.Vary2Edit.Position = [7 165 217 22];

            % Create Vary2ListBox
            app.Vary2ListBox = uilistbox(app.Vary2Tab);
            app.Vary2ListBox.Items = {};
            app.Vary2ListBox.ValueChangedFcn = createCallbackFcn(app, @Vary2List, true);
            app.Vary2ListBox.Position = [7 8 271 150];
            app.Vary2ListBox.Value = {};

            % Create Vary2RemoveButton
            app.Vary2RemoveButton = uibutton(app.Vary2Tab, 'push');
            app.Vary2RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @Vary2Remove, true);
            app.Vary2RemoveButton.FontSize = 8;
            app.Vary2RemoveButton.Position = [281 128 44 18];
            app.Vary2RemoveButton.Text = 'Remove';

            % Create Vary3Tab
            app.Vary3Tab = uitab(app.VaryTabGroup);
            app.Vary3Tab.AutoResizeChildren = 'off';
            app.Vary3Tab.Title = 'Vary3';

            % Create Vary3Switch
            app.Vary3Switch = uiswitch(app.Vary3Tab, 'slider');
            app.Vary3Switch.Position = [253 166 45 20];
            app.Vary3Switch.Value = 'On';

            % Create Vary3Edit
            app.Vary3Edit = uieditfield(app.Vary3Tab, 'text');
            app.Vary3Edit.ValueChangedFcn = createCallbackFcn(app, @EditVary3, true);
            app.Vary3Edit.Position = [7 165 217 22];

            % Create Vary3ListBox
            app.Vary3ListBox = uilistbox(app.Vary3Tab);
            app.Vary3ListBox.Items = {};
            app.Vary3ListBox.ValueChangedFcn = createCallbackFcn(app, @Vary3List, true);
            app.Vary3ListBox.Position = [7 8 271 150];
            app.Vary3ListBox.Value = {};

            % Create Vary3RemoveButton
            app.Vary3RemoveButton = uibutton(app.Vary3Tab, 'push');
            app.Vary3RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @Vary3Remove, true);
            app.Vary3RemoveButton.FontSize = 8;
            app.Vary3RemoveButton.Position = [281 128 44 18];
            app.Vary3RemoveButton.Text = 'Remove';

            % Create Vary4Tab
            app.Vary4Tab = uitab(app.VaryTabGroup);
            app.Vary4Tab.AutoResizeChildren = 'off';
            app.Vary4Tab.Title = 'Vary4';

            % Create Vary4Switch
            app.Vary4Switch = uiswitch(app.Vary4Tab, 'slider');
            app.Vary4Switch.Position = [253 166 45 20];
            app.Vary4Switch.Value = 'On';

            % Create Vary4Edit
            app.Vary4Edit = uieditfield(app.Vary4Tab, 'text');
            app.Vary4Edit.ValueChangedFcn = createCallbackFcn(app, @EditVary4, true);
            app.Vary4Edit.Position = [7 165 217 22];

            % Create Vary4ListBox
            app.Vary4ListBox = uilistbox(app.Vary4Tab);
            app.Vary4ListBox.Items = {};
            app.Vary4ListBox.ValueChangedFcn = createCallbackFcn(app, @Vary4List, true);
            app.Vary4ListBox.Position = [7 8 271 150];
            app.Vary4ListBox.Value = {};

            % Create Vary4RemoveButton
            app.Vary4RemoveButton = uibutton(app.Vary4Tab, 'push');
            app.Vary4RemoveButton.ButtonPushedFcn = createCallbackFcn(app, @Vary4Remove, true);
            app.Vary4RemoveButton.FontSize = 8;
            app.Vary4RemoveButton.Position = [281 128 44 18];
            app.Vary4RemoveButton.Text = 'Remove';

            % Create FitOptTab
            app.FitOptTab = uitab(app.VaryTabGroup);
            app.FitOptTab.AutoResizeChildren = 'off';
            app.FitOptTab.Title = 'FitOpt';

            % Create FitOptEdit
            app.FitOptEdit = uieditfield(app.FitOptTab, 'text');
            app.FitOptEdit.ValueChangedFcn = createCallbackFcn(app, @EditFitOpt, true);
            app.FitOptEdit.Position = [178 163 142 22];

            % Create MethodDropDownLabel
            app.MethodDropDownLabel = uilabel(app.FitOptTab);
            app.MethodDropDownLabel.HorizontalAlignment = 'right';
            app.MethodDropDownLabel.Position = [7 164 47 22];
            app.MethodDropDownLabel.Text = 'Method';

            % Create FitMethodDropDown
            app.FitMethodDropDown = uidropdown(app.FitOptTab);
            app.FitMethodDropDown.Items = {'simplex', 'levmar', 'montecarlo', 'genetic', 'grid', 'swarm'};
            app.FitMethodDropDown.Position = [69 164 89 22];
            app.FitMethodDropDown.Value = 'simplex';

            % Create TargetDropDownLabel
            app.TargetDropDownLabel = uilabel(app.FitOptTab);
            app.TargetDropDownLabel.HorizontalAlignment = 'right';
            app.TargetDropDownLabel.Position = [16 134 38 22];
            app.TargetDropDownLabel.Text = 'Target';

            % Create TargetDropDown
            app.TargetDropDown = uidropdown(app.FitOptTab);
            app.TargetDropDown.Items = {'as is', 'integral', 'dbl integral', 'derivative', 'FFT', ''};
            app.TargetDropDown.Position = [69 134 89 22];
            app.TargetDropDown.Value = 'as is';

            % Create ScalingDropDownLabel
            app.ScalingDropDownLabel = uilabel(app.FitOptTab);
            app.ScalingDropDownLabel.HorizontalAlignment = 'right';
            app.ScalingDropDownLabel.Position = [9 103 45 22];
            app.ScalingDropDownLabel.Text = 'Scaling';

            % Create FitScalingDropDown
            app.FitScalingDropDown = uidropdown(app.FitOptTab);
            app.FitScalingDropDown.Items = {'minmax', 'maxabs', 'lsq', 'lsq0', 'lsq1', 'lsq2', 'none', ''};
            app.FitScalingDropDown.Position = [69 103 89 22];
            app.FitScalingDropDown.Value = 'minmax';

            % Create DDataButtonGroup
            app.DDataButtonGroup = uibuttongroup(app.FitOptTab);
            app.DDataButtonGroup.AutoResizeChildren = 'off';
            app.DDataButtonGroup.Title = '2D Data';
            app.DDataButtonGroup.Position = [26 9 114 85];

            % Create AllParallelButton
            app.AllParallelButton = uiradiobutton(app.DDataButtonGroup);
            app.AllParallelButton.Text = 'All Parallel';
            app.AllParallelButton.Position = [6 23 78 22];

            % Create AllSequentialButton
            app.AllSequentialButton = uiradiobutton(app.DDataButtonGroup);
            app.AllSequentialButton.Text = 'All Sequential';
            app.AllSequentialButton.Position = [6 3 95 22];

            % Create OnlyCurrentButton
            app.OnlyCurrentButton = uiradiobutton(app.DDataButtonGroup);
            app.OnlyCurrentButton.Text = 'Only Current';
            app.OnlyCurrentButton.Position = [6 42 90 22];
            app.OnlyCurrentButton.Value = true;

            % Create FitOptListBox
            app.FitOptListBox = uilistbox(app.FitOptTab);
            app.FitOptListBox.Items = {};
            app.FitOptListBox.ValueChangedFcn = createCallbackFcn(app, @FitOptList, true);
            app.FitOptListBox.Position = [178 7 142 149];
            app.FitOptListBox.Value = {};

            % Create SimulateButton
            app.SimulateButton = uibutton(app.cwEPRapp, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @Simulate, true);
            app.SimulateButton.FontSize = 14;
            app.SimulateButton.Position = [933 38 84 28];
            app.SimulateButton.Text = 'Simulate';

            % Create Spinner
            app.Spinner = uispinner(app.cwEPRapp);
            app.Spinner.Limits = [1 Inf];
            app.Spinner.ValueChangedFcn = createCallbackFcn(app, @spinner, true);
            app.Spinner.Visible = 'off';
            app.Spinner.Position = [783 305 56 22];
            app.Spinner.Value = 1;

            % Create SliceLabel
            app.SliceLabel = uilabel(app.cwEPRapp);
            app.SliceLabel.Visible = 'off';
            app.SliceLabel.Position = [796 326 31 22];
            app.SliceLabel.Text = 'Slice';

            % Create DataFilenameLabel
            app.DataFilenameLabel = uilabel(app.cwEPRapp);
            app.DataFilenameLabel.Position = [113 682 630 22];
            app.DataFilenameLabel.Text = 'DataFilename';

            % Create gvaluesButton
            app.gvaluesButton = uibutton(app.cwEPRapp, 'state');
            app.gvaluesButton.ValueChangedFcn = createCallbackFcn(app, @gvalues, true);
            app.gvaluesButton.Text = 'g values';
            app.gvaluesButton.Position = [11 236 66 22];

            % Create SimTimeLabel
            app.SimTimeLabel = uilabel(app.cwEPRapp);
            app.SimTimeLabel.HorizontalAlignment = 'right';
            app.SimTimeLabel.FontSize = 10;
            app.SimTimeLabel.Visible = 'off';
            app.SimTimeLabel.Position = [908 14 107 16];
            app.SimTimeLabel.Text = 'Sim Time =';

            % Create FittingSliceLabel
            app.FittingSliceLabel = uilabel(app.cwEPRapp);
            app.FittingSliceLabel.Visible = 'off';
            app.FittingSliceLabel.Position = [504 236 323 22];
            app.FittingSliceLabel.Text = 'Fitting Slice ';

            % Create StopFittingButton
            app.StopFittingButton = uibutton(app.cwEPRapp, 'state');
            app.StopFittingButton.ValueChangedFcn = createCallbackFcn(app, @StopFitting, true);
            app.StopFittingButton.Visible = 'off';
            app.StopFittingButton.Text = 'Stop Fitting';
            app.StopFittingButton.Position = [738 236 89 22];

            % Create esfitButton
            app.esfitButton = uibutton(app.cwEPRapp, 'push');
            app.esfitButton.ButtonPushedFcn = createCallbackFcn(app, @esfit, true);
            app.esfitButton.FontSize = 14;
            app.esfitButton.Position = [844 38 83 28];
            app.esfitButton.Text = 'esfit';

            % Create FunctionDropDown
            app.FunctionDropDown = uidropdown(app.cwEPRapp);
            app.FunctionDropDown.Items = {'garlic', 'chili', 'pepper'};
            app.FunctionDropDown.ValueChangedFcn = createCallbackFcn(app, @FunctionSelection, true);
            app.FunctionDropDown.Position = [1026 40 79 22];
            app.FunctionDropDown.Value = 'garlic';

            % Create SimMethodDropDown
            app.SimMethodDropDown = uidropdown(app.cwEPRapp);
            app.SimMethodDropDown.Items = {'exact', 'perturb', 'perturb1', 'perturb2', 'perturb3', 'perturb4', 'perturb5', ''};
            app.SimMethodDropDown.Position = [1026 14 79 22];
            app.SimMethodDropDown.Value = 'perturb2';

            % Create FitGUICheckBox
            app.FitGUICheckBox = uicheckbox(app.cwEPRapp);
            app.FitGUICheckBox.Text = 'Fit GUI';
            app.FitGUICheckBox.Position = [849 11 60 22];
            app.FitGUICheckBox.Value = true;

            % Create RescaleDropDown
            app.RescaleDropDown = uidropdown(app.cwEPRapp);
            app.RescaleDropDown.Items = {'minmax', 'maxabs', 'lsq', 'lsq0', 'lsq1', 'lsq2', 'none'};
            app.RescaleDropDown.ValueChangedFcn = createCallbackFcn(app, @showDatas, true);
            app.RescaleDropDown.Position = [1112 14 79 22];
            app.RescaleDropDown.Value = 'minmax';

            % Create RescaleLabel
            app.RescaleLabel = uilabel(app.cwEPRapp);
            app.RescaleLabel.Position = [1128 35 48 22];
            app.RescaleLabel.Text = 'Rescale';

            % Create v33Label
            app.v33Label = uilabel(app.cwEPRapp);
            app.v33Label.Position = [1170 696 28 22];
            app.v33Label.Text = 'v3.3';

            % Create RealCheckBox
            app.RealCheckBox = uicheckbox(app.cwEPRapp);
            app.RealCheckBox.ValueChangedFcn = createCallbackFcn(app, @RealCheck, true);
            app.RealCheckBox.Visible = 'off';
            app.RealCheckBox.Text = 'Real';
            app.RealCheckBox.Position = [87 236 46 22];

            % Create ImaginaryCheckBox
            app.ImaginaryCheckBox = uicheckbox(app.cwEPRapp);
            app.ImaginaryCheckBox.ValueChangedFcn = createCallbackFcn(app, @ImagCheck, true);
            app.ImaginaryCheckBox.Visible = 'off';
            app.ImaginaryCheckBox.Text = 'Imaginary';
            app.ImaginaryCheckBox.Position = [135 236 74 22];
        end
    end

    methods (Access = public)

        % Construct app
        function app = cwEPR

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.cwEPRapp)

            % Execute the startup function
            runStartupFcn(app, @Opening)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.cwEPRapp)
        end
    end
end
