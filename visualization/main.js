// You can only use import and export statements inside modules; not regular scripts.
// In some module systems, you can omit the file extension and the dot (e.g. '/modules/square').
// This doesn't work in native JavaScript modules.
import './chart.js'

// https://stackoverflow.com/a/39855320/6674256
import regions from '../data/sample_of_regions.js'

new Vue({
  el: '#app',
  vuetify: new Vuetify(),
  data () {
    return {
      regions: regions,
      names: Object.keys(regions),
      name: Object.keys(regions)[0],
      colors: ['red', 'green', 'blue', 'magenta', 'cyan'],
      desserts: [
        {
          name: 'Frozen Yogurt',
          calories: 159,
        },
        {
          name: 'Ice cream sandwich',
          calories: 237,
        },
        {
          name: 'Eclair',
          calories: 262,
        },
        {
          name: 'Cupcake',
          calories: 305,
        },
        {
          name: 'Gingerbread',
          calories: 356,
        },
        {
          name: 'Jelly bean',
          calories: 375,
        },
        {
          name: 'Lollipop',
          calories: 392,
        },
        {
          name: 'Honeycomb',
          calories: 408,
        },
        {
          name: 'Donut',
          calories: 452,
        },
        {
          name: 'KitKat',
          calories: 518,
        },
      ],
    }
  },
  // use regular functions, not arrow functions, when invoking "this" : https://michaelnthiessen.com/this-is-undefined/
  computed: {
    region () {
      return this.regions[this.name]
    },
    number_samples () {
      return this.region['normalized read depth'].length
    },
    // plotly only supports a limited number of html tags, e.g., "<br>"
    hover_text () {
      return this.region['inferred copy number'].map((e, i) => {
        const first_line = 'inferred copy number: ' + e + '<br>'
        const second_line = 'possible copy numbers: ' + this.region['possible copy numbers'] + '<br>'
        const third_line = ('probabilities of possible copy numbers: ' +
          this.region['probabilities of possible copy numbers'][i].map(value => value.toFixed(2)))
        return first_line + second_line + third_line
      })
    },
    data () {
      return {
        traces: [
          {
            x: this.region['normalized read depth'],
            y: Array.from({ length: this.number_samples }, () => Math.random()),
            mode: 'markers',
            type: 'scatter',
            marker: {
              size: 12,
              color: this.region['inferred copy number'].map(value => this.colors[value - 1])
            },
            text: this.hover_text,
            hoverinfo: 'text'
          }
        ],
        layout: {
          font: {
            family: 'Google Sans, sans-serif',
          },
          xaxis: {
            title: 'normalized read depth'
          },
          // https://codepen.io/plotly/pen/KpLVzv ...
          yaxis: {
            showgrid: false,
            zeroline: false,
            showticklabels: false
          },
          hovermode: 'closest',
          hoverlabel: {
            bgcolor: '#fafafa'
          }
        }
      }
    }
  },
  methods: {
    print (item) {
      console.log('item is ', item.name)
    },
  },
  mounted () {
    console.log(this.regions)
  },
})
