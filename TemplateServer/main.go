package main

import (
	"archive/tar"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"sort"
	"strings"
	"time"

	"github.com/goji/httpauth"
	"github.com/gorilla/handlers"
	"github.com/gorilla/mux"
	"github.com/rs/cors"
)

func ParseConfigName(args []string) (string, []string) {
	resArgs := make([]string, 0)
	file := ""
	for i := 0; i < len(args); i++ {
		if args[i] == "-config" {
			file = args[i+1]
			i++
			continue
		}

		resArgs = append(resArgs, args[i])
	}

	return file, resArgs
}

func dbpaths(path string) (string, string) {
	return path + ".ffdata", path + ".ffindex"
}

func unique(slice []string) []string {
	check := make(map[string]bool, 0)
	unique := make([]string, 0)
	for _, e := range slice {
		if len(e) != 0 && !check[e] {
			unique = append(unique, e)
			check[e] = true
		}
	}
	return unique
}

func AddTarEntry(tw *tar.Writer, name string, s string, now time.Time) error {
	header := new(tar.Header)
	header.Name = name
	header.Size = int64(len(s))
	header.Mode = int64(0666)
	header.ModTime = now
	if err := tw.WriteHeader(header); err != nil {
		return err
	}
	if _, err := tw.Write([]byte(s)); err != nil {
		return err
	}
	return nil
}

func GatherResults(w io.Writer, templates []string, a3m Reader, hhm Reader, cif Reader) (err error) {
	gw := gzip.NewWriter(w)
	tw := tar.NewWriter(gw)

	uniques := make([]string, len(templates))
	uniquesWithoutChains := make([]string, len(templates))
	for i := 0; i < len(templates); i++ {
		key := templates[i]
		slice := strings.Split(key, "_")
		if len(slice) != 2 {
			continue
		}
		uniques = append(uniques, key)
		uniquesWithoutChains = append(uniquesWithoutChains, slice[0])
	}
	uniques = unique(uniques)
	uniquesWithoutChains = unique(uniquesWithoutChains)

	sort.Strings(uniques)
	sort.Strings(uniquesWithoutChains)

	a3mOffset := 0
	a3mData := strings.Builder{}
	a3mIndex := strings.Builder{}

	hhmOffset := 0
	hhmData := strings.Builder{}
	hhmIndex := strings.Builder{}
	for i := 0; i < len(uniques); i++ {
		a3mid, ok := a3m.Id(uniques[i])
		if ok == false {
			continue
		}
		hhmid, ok := hhm.Id(uniques[i])
		if ok == false {
			continue
		}
		entry := a3m.Data(a3mid)
		entryLen := len(entry) + 1
		a3mData.WriteString(entry)
		a3mData.WriteRune(rune(0))
		a3mIndex.WriteString(fmt.Sprintf("%s\t%d\t%d\n", uniques[i], a3mOffset, entryLen))
		a3mOffset += entryLen

		entry = hhm.Data(hhmid)
		entryLen = len(entry) + 1
		hhmData.WriteString(entry)
		hhmData.WriteRune(rune(0))
		hhmIndex.WriteString(fmt.Sprintf("%s\t%d\t%d\n", uniques[i], hhmOffset, entryLen))
		hhmOffset += entryLen
	}

	now := time.Now()
	if err := AddTarEntry(tw, "pdb70_a3m.ffdata", a3mData.String(), now); err != nil {
		return err
	}
	if err := AddTarEntry(tw, "pdb70_a3m.ffindex", a3mIndex.String(), now); err != nil {
		return err
	}
	if err := AddTarEntry(tw, "pdb70_hhm.ffdata", hhmData.String(), now); err != nil {
		return err
	}
	if err := AddTarEntry(tw, "pdb70_hhm.ffindex", hhmIndex.String(), now); err != nil {
		return err
	}

	for i := 0; i < len(uniquesWithoutChains); i++ {
		pdbacc := uniquesWithoutChains[i]
		cifid, ok := cif.Id(pdbacc)
		if ok == false {
			continue
		}

		if err := AddTarEntry(tw, strings.ToLower(pdbacc)+".cif", cif.Data(cifid), now); err != nil {
			return err
		}
	}

	if err := tw.Close(); err != nil {
		return err
	}

	if err := gw.Close(); err != nil {
		return err
	}

	return nil
}

func main() {
	configFile, args := ParseConfigName(os.Args[1:])

	var config ConfigRoot
	var err error
	if len(configFile) > 0 {
		config, err = ReadConfigFromFile(configFile)

	} else {
		config, err = DefaultConfig()
	}
	if err != nil {
		panic(err)
	}

	err = config.ReadParameters(args)
	if err != nil {
		panic(err)
	}

	baseRouter := mux.NewRouter()
	var r *mux.Router
	if len(config.Server.PathPrefix) > 0 {
		r = baseRouter.PathPrefix(config.Server.PathPrefix).Subrouter()
	} else {
		r = baseRouter
	}

	a3mreader := Reader{}
	a3mreader.Make(dbpaths(config.Paths.DatabasePrefix + "_a3m"))

	hhmreader := Reader{}
	hhmreader.Make(dbpaths(config.Paths.DatabasePrefix + "_hhm"))

	cifreader := Reader{}
	cifreader.Make(dbpaths(config.Paths.DatabasePrefix + "_cif"))

	r.HandleFunc("/template/{list}", func(w http.ResponseWriter, req *http.Request) {
		templates := strings.Split(mux.Vars(req)["list"], ",")
		w.Header().Set("Content-Disposition", "attachment; filename=\"templates.tar.gz\"")
		w.Header().Set("Content-Type", "application/octet-stream")
		if err != GatherResults(w, templates, a3mreader, hhmreader, cifreader) {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
	}).Methods("GET")

	h := http.Handler(r)
	if config.Server.Auth != nil {
		h = httpauth.SimpleBasicAuth(config.Server.Auth.Username, config.Server.Auth.Password)(h)
	}
	if config.Verbose {
		h = handlers.LoggingHandler(os.Stdout, h)
	}
	if config.Server.CORS {
		c := cors.AllowAll()
		h = c.Handler(h)
	}

	srv := &http.Server{
		Handler: h,
		Addr:    config.Server.Address,
	}

	log.Println("MMseqs2 Webserver")
	log.Fatal(srv.ListenAndServe())
}
