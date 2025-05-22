import pandas as pd
import altair as alt


class WACIAttribtion:
    portfolio_t1_df = pd.DataFrame()
    portfolio_t2_df = None
    esg_data_t1_df = pd.DataFrame()
    esg_data_t2_df = None
    pf_esg_t1_df = pd.DataFrame
    pf_esg_t2_df = None
    identifier = ""
    nmv_name = ""
    carbon_intensity_name = ""
    config_def = []
    waci_t1 = None
    waci_t2 = None
    attribution = None

    def __init__(self, config_definition):
        #Disable chain assignment warning
        pd.options.mode.chained_assignment = None

        #Initialize class parameters
        self.config_def = config_definition
        self.carbon_intensity_name = self.config_def['carbon_intensity']
        self.pf_identifier = self.config_def['pf_identifier']
        self.esg_identifier = self.config_def['esg_identifier']
        self.nmv_name = self.config_def['nmv_name']

        #Read portfolio and ESG data as time t
        self.portfolio_t1_df = pd.read_csv(self.config_def['pf_t1'])
        self.esg_data_t1_df = pd.read_csv(self.config_def['esg_t1'])
        #Merge portfolio and esg data on identifier
        self.pf_esg_t1_df = pd.merge(self.portfolio_t1_df, self.esg_data_t1_df, left_on=self.pf_identifier,
                                     right_on=self.esg_identifier,
                                     how='left')
        # Read portfolio and ESG data as time t+1
        if 'pf_t2' in config_definition.keys():
            self.portfolio_t2_df = pd.read_csv(self.config_def['pf_t2'])
        else:
            self.portfolio_t2_df = None
            return

        if 'esg_t2' in config_definition.keys():
            self.esg_data_t2_df = pd.read_csv(self.config_def['esg_t2'])
        else:
            self.esg_data_t2_df = None
            return
        # Merge portfolio and esg data on identifier
        self.pf_esg_t2_df = pd.merge(self.portfolio_t2_df, self.esg_data_t2_df, left_on=self.pf_identifier,
                                     right_on=self.esg_identifier,
                                     how='left')
        return

    def calculate_waci(self, in_scope_col):
        """
            Calculates Weighted Average Carbon Intensity (WACI) of the 2 portfolios at time t and time t-1

            This function takes the loaded portfolios during initialization and compute the WACI of the 2 portfolios.
            It uses the in_scope_col parameter and filter out securities that are in-scope before calculating.

            Args:
                in_scope_col (string):  The column header of the columns that indicate in-scope assets. A value of 1
                                        indicates asset is sin-scope.

            Returns:
                string,string: WACI at t and at t-1

            Raises:
                None
            """

        #Determine in scope asset using in_scope_col value
        pf_t1 = self.pf_esg_t1_df[self.pf_esg_t1_df[in_scope_col] == 1]
        pf_t2 = self.pf_esg_t2_df[self.pf_esg_t1_df[in_scope_col] == 1]

        #Calculates WACI at t and at t+1 if the data is loaded during initialization
        self.waci_t1 = (pf_t1[self.nmv_name] * pf_t1[self.carbon_intensity_name]
                        /pf_t1.dropna(subset=[self.carbon_intensity_name])[
                            self.nmv_name].sum()).sum()
        if self.pf_esg_t2_df is not None:
            self.waci_t2 = (pf_t2[self.nmv_name] * pf_t2[self.carbon_intensity_name]
                            /pf_t2.dropna(subset=[self.carbon_intensity_name])[self.nmv_name].sum()).sum()

        return self.waci_t1, self.waci_t2

    def calculate_attribution(self, in_scope_col, exclusion_col=""):
        """
            Calculates the WACI attribution between the 2 portfolio at time t and time t-1

            This function takes the loaded portfolios during initialization and compute the WACI attribution between
            time t and time t-1. It uses the in_scope_col parameter and filter out securities that are in-scope before
            calculating. It is following the NZAOA calculation methodology and spits out attribution to:
            Coverage, New Investments, Divestments, Weights, Emissions, Revenue and Model.
            Model attribution happens when Carbon Intensity <> Carbon Emissions/Sales. This is when data provider
            models the intensity rather than using this standard formula.

            Args:
                in_scope_col (string):  The column header of the columns that indicate in-scope assets. A value of 1
                                        indicates in-scope.
                exclusion_col (string): Intended to support asset management specific exclusions. Not used.
            Returns:
                dict: WACI attributions

            Raises:
                None
            """
        attribution_dic = {}

        #Calculates portfolio weights of in-scope assets at t
        pf_t1 = self.pf_esg_t1_df[self.pf_esg_t1_df[in_scope_col] == 1]
        pf_t2 = self.pf_esg_t2_df[self.pf_esg_t1_df[in_scope_col] == 1]
        pf_t1['weight'] = pf_t1[self.nmv_name] / pf_t1.dropna(subset=[self.carbon_intensity_name])[self.nmv_name].sum()

        #Calculates portfolio weights of in-scope assets at t+1
        pf_t2['weight'] = pf_t2[self.nmv_name] / pf_t2.dropna(subset=[self.carbon_intensity_name])[self.nmv_name].sum()

        #Calculates attributiononly when there is value WACI at t and t+1
        if self.waci_t1 != None and self.waci_t2 != None:
            # Calculate coverage attribution
            #pf_t1_no_esg_data_df = pf_t1[pf_t1[self.carbon_intensity_name].isna()]
            #pf_t2_has_esg_data_df = pf_t2.dropna(subset=[self.carbon_intensity_name])
            #new_coverage = pf_t2_has_esg_data_df[
            #    pf_t2_has_esg_data_df[self.pf_identifier].isin(pf_t1_no_esg_data_df[self.pf_identifier])]

            pf_t1_t2_df = pd.merge(pf_t1, pf_t2, on=self.pf_identifier, how='outer')

            #Calculates attribution due to new Investment
            new_investment_df = pf_t1_t2_df[pf_t1_t2_df[self.nmv_name + "_x"].isna()]
            new_investment_df['waci_attri'] = new_investment_df['weight_y'] * (
                    new_investment_df[self.carbon_intensity_name + "_y"] - self.waci_t1)
            attri_new = new_investment_df['waci_attri'].sum()
            #Calculates attribution due to divestment
            divestment_df = pf_t1_t2_df[pf_t1_t2_df[self.nmv_name + "_y"].isna()]
            divestment_df['waci_attri'] = divestment_df['weight_x'] * (
                    divestment_df[self.carbon_intensity_name + "_x"] - self.waci_t1)
            attri_div = divestment_df['waci_attri'].sum()

            #Calculates attribution for held security - Weight, Emissions, Sales, Model
            intersect_df = pf_t1_t2_df.dropna(subset=[self.nmv_name + "_x", self.nmv_name + "_y"])

            #intersect - extract security without standard intensity calculation
            #check that it is a custom model at t1 and has a valid data at t2. It is a standard WACI
            #if it use Carbon Emissions/Sales, otherwise it is non-standard.
            intersect_df['inten_check'] = intersect_df[self.config_def['carbon_emission'] + "_x"] / intersect_df[
                self.config_def['sales'] + "_x"]
            custom_model_df = intersect_df[
                (intersect_df['inten_check'] != intersect_df[self.carbon_intensity_name + "_x"]) &
                (intersect_df[self.carbon_intensity_name + "_y"] > 0)]
            # drop all other rows that has missing data
            custom_model_df = custom_model_df.dropna(subset=['inten_check'])
            custom_model_df['waci_attri'] = custom_model_df['weight_y'] * (
                    custom_model_df[self.carbon_intensity_name + "_y"] - custom_model_df[
                self.carbon_intensity_name + "_x"])
            attri_custom = custom_model_df['waci_attri'].sum()

            #Calculates attribution due to coverage
            new_coverage_df = intersect_df[
                (intersect_df[self.nmv_name + "_x"] > 0) & (intersect_df[self.carbon_intensity_name + "_x"].isna())]
            new_coverage_df['waci_attri'] = new_coverage_df['weight_y'] * (
                    new_coverage_df[self.carbon_intensity_name + "_y"] - self.waci_t1)
            attri_new_coverage = new_coverage_df['waci_attri'].sum()

            reduce_coverage_df = intersect_df[
                (intersect_df[self.nmv_name + "_y"] > 0) & (intersect_df[self.carbon_intensity_name + "_y"].isna())]
            reduce_coverage_df['waci_attri'] = reduce_coverage_df['weight_x'] * (
                    reduce_coverage_df[self.carbon_intensity_name + "_x"] - self.waci_t1)
            attri_reduce_coverage = reduce_coverage_df['waci_attri'].sum()

            #Calculates weights, revenue, emissions attributions
            intersect_df = intersect_df[~intersect_df[self.pf_identifier].isin(new_coverage_df[self.pf_identifier])]
            intersect_df = intersect_df[~intersect_df[self.pf_identifier].isin(reduce_coverage_df[self.pf_identifier])]
            intersect_df['weight_attri'] = (intersect_df[self.carbon_intensity_name + "_x"] - self.waci_t1) * \
                                           (intersect_df['weight_y'] - intersect_df['weight_x'])
            intersect_df['rev_attri'] = intersect_df['weight_y'] * \
                                        (intersect_df[self.config_def['carbon_emission'] + "_y"] + intersect_df[
                                            self.config_def['carbon_emission'] + "_x"]) / 2 * \
                                        (1 / intersect_df[self.config_def['sales'] + "_y"] - 1 / intersect_df[
                                            self.config_def['sales'] + "_x"])
            intersect_df['carbon_attri'] = intersect_df['weight_y'] * \
                                           (intersect_df[self.config_def['carbon_emission'] + "_y"] - intersect_df[
                                               self.config_def['carbon_emission'] + "_x"]) * \
                                           (1 / intersect_df[self.config_def['sales'] + "_y"] + 1 / intersect_df[
                                               self.config_def['sales'] + "_x"]) / 2
            attri_rev = intersect_df['rev_attri'].sum()
            attri_carbon = intersect_df['carbon_attri'].sum()
            attri_weight = intersect_df['weight_attri'].sum()

            #Assign attribution to return dictionary
            attribution_dic['Weight'] = attri_weight
            attribution_dic['Coverage'] = attri_new_coverage - attri_reduce_coverage
            attribution_dic['Model'] = attri_custom
            attribution_dic['Investment'] = attri_new
            attribution_dic['Divestment'] = -attri_div
            attribution_dic['Emissions'] = attri_carbon
            attribution_dic['Revenue'] = attri_rev

            #check for error, sum of all attribution and old WACI should result in new WACI value
            error = self.waci_t1 + attri_weight + attri_new_coverage - attri_reduce_coverage + attri_custom \
                    + attri_new - attri_div + attri_carbon + attri_rev - self.waci_t2
            print(error)
            self.attribution = attribution_dic
            return attribution_dic

    def plot_attribution(self,
                         order=['Coverage', 'Model', 'Investment', 'Divestment', 'Weight', 'Emissions', 'Revenue']):
        """
            Plot the WACI attribution as a cumulative chart

            Args:
                order (list): The order which the attributions will appear in the chart from left to right
            Returns:
                chart (Altair chart): WACI attributions cumulative chart

            Raises:
                None
        """

        # create dataframe for use in chart
        chart_df = pd.DataFrame([self.attribution])[order].T.reset_index()
        chart_df.columns = ['Attribution', 'Value']

        #calculates the range of the bar chart
        chart_df['High'] = chart_df.cumsum()['Value'] + self.waci_t1
        chart_df['Low'] = chart_df['High'] - chart_df['Value']

        #generate the direction value on whether it is a increase or decrease from previous bar
        chart_df['Direction'] = ['Decrease' if row['Low'] > row['High'] else 'Increase' for idx, row in
                                 chart_df[['High', 'Low']].iterrows()]
        min_waci = chart_df[['High', 'Low']].min().min()
        max_waci = chart_df[['High', 'Low']].max().max()

        #add in the WACI values before and after attribution and assign custom direction values
        chart_df = pd.concat(
            [pd.DataFrame([['Old WACI', self.waci_t1, min_waci * 0.8, self.waci_t1, 'Old WACI']],
                          columns=['Attribution', 'High', 'Low', 'Value', 'Direction']), chart_df,
             pd.DataFrame([['New WACI', self.waci_t2, min_waci * 0.8, self.waci_t2, 'New WACI']],
                          columns=['Attribution', 'High', 'Low', 'Value', 'Direction'])
             ])
        #Generate the bar chart with custom colors
        chart = alt.Chart(chart_df).mark_bar().encode(
            x=alt.X('Attribution:O', title='Attribution Factors',
                    sort=['Old WACI'] + order + ['New WACI']),
            y=alt.Y('Low:Q', title='WACI', scale=alt.Scale(domain=[min_waci * 0.8, max_waci])),
            y2='High:Q',
            color=alt.Color('Direction:O',
                            scale=alt.Scale(domain=['Increase', 'Decrease', 'Old WACI',
                                                    'New WACI'],
                                            range=['red', 'green', '#116EA1', '#49A6D9'])),
            tooltip=['Low', 'High']).properties(title='WACI Attribution').configure(background='#D9E9F0')

        return chart
