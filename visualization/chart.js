// https://medium.com/@bobthomas295/vuejs-and-plotly-js-reactive-charts-da9b3b59f2dc
Vue.component('chart', {
  props: ['data'],
  template: '<div ref="chart"></div> ',
  mounted () {
    Plotly.plot(this.$refs.chart, this.data.traces, this.data.layout, {responsive: true})
  },
  watch: {
    data: {
      handler: function () {
        Plotly.react(
          this.$refs.chart,
          this.data.traces,
          this.data.layout
        )
      },
      deep: true
    }
  }
})

