function [post_z,pre_z]=stan_analyze_nervecut_audio()
%
%
%

% get SAP features, compare pre/post for best match motif
%


save_dir='features';

% get options

[options,dirs]=stan_preflight;

key=stan_read_nervecut_audio;

feature_dir=fullfile(dirs.agg_dir,dirs.nervecut_audio_dir,save_dir);
listing=dir(fullfile(feature_dir,'*.mat'));
listing={listing(:).name};
post_z=struct();

for i=1:length(listing)

	listing{i}
	load(fullfile(feature_dir,listing{i}),'pre_features','post_features');

	feature_names=fieldnames(pre_features);
	feature_names(strcmp(feature_names,'spec_deriv'))=[];

	pre_amp=mean(zscore(pre_features.amp)');
	post_amp=mean(zscore(post_features.amp)');

	pre_amp_idx=pre_amp>=options.audio_amp_thresh;
	post_amp_idx=post_amp>=options.audio_amp_thresh;

	for j=1:length(feature_names)

		mu1=mean(pre_features.(feature_names{j})(pre_amp_idx,:));
		var1=var(pre_features.(feature_names{j})(pre_amp_idx,:));

		% convert post to z

		mu2=mean(post_features.(feature_names{j})(post_amp_idx,:));
		var2=var(post_features.(feature_names{j})(post_amp_idx,:));

		post_z.ranksum_mu(i).(feature_names{j})=ranksum(mu1,mu2);
		post_z.ranksum_var(i).(feature_names{j})=ranksum(var1,var2);

		post_z.mu(i).(feature_names{j})=(mu2(1:100)-mean(mu1))/std(mu1);
		post_z.var(i).(feature_names{j})=(var2(1:100)-mean(var1))/std(var1);

		pre_z.raw_mu(i).(feature_names{j})=mu1(1:100);
		pre_z.raw_var(i).(feature_names{j})=var1(1:100);
		post_z.raw_mu(i).(feature_names{j})=mu2(1:100);
		post_z.raw_var(i).(feature_names{j})=var2(1:100);

	end

end



%%%% print out table with ranksum statistics

fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'fig3_audionervecut.txt'),'w+');

for i=1:length(feature_names)

	[p_mu,h,stats_mu]=ranksum(cat(2,pre_z(:).raw_mu.(feature_names{i})),cat(2,post_z(:).raw_mu.(feature_names{i})));
	[p_var,h,stats_var]=ranksum(cat(2,pre_z(:).raw_var.(feature_names{i})),cat(2,post_z(:).raw_var.(feature_names{i})));

	fprintf(fid,'%-20s\tp=%-2.3e\tz=%g\t%-20s\tp=%-2.3e\tz=%g\n',...
		[ feature_names{i} '(mean)' ] ,p_mu,stats_mu.zval,...
		[ feature_names{i} '(var)' ],p_var,stats_var.zval);
end

fprintf(fid,'N(points): %i',length(pre_z.raw_mu(1).(feature_names{1})));
fclose(fid);
