// https://medium.com/@bobthomas295/vuejs-and-plotly-js-reactive-charts-da9b3b59f2dc
Vue.component('chart', {
  props: ['data'],
  template: '<div ref="chart" class="elevation-10"></div>',
  data () {
    return {
      tippy: null,
      tippy_timer: null,
      wait: 1000
    }
  },
  mounted () {
    Plotly.plot(this.$refs.chart, this.data.traces, this.data.layout, {responsive: true})
    this.tippy_timer = setTimeout(this.show_tippy, this.wait)
  },
  methods: {
    show_tippy () {
      this.tippy = tippy(this.$refs.chart, {
        content: 'Hover over samples in scatter plot to see confidence in inferred copy number',
        arrow: true,
        theme: 'google',
        trigger: 'manual',
        hideOnClick: true
      })
      this.tippy.show()
    }
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
  },
  beforeDestroy () {
    clearTimeout(this.tippy_timer)
    this.tippy.destroy(true)
  }
})

