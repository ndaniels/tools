export GOBIN=./bin

all:
	go install ./...

clean:
	rm -f bin/*

fmt:
	gofmt -w *.go */*.go
	colcheck *.go */*.go

push:
	git push origin master
	git push github master
	git push tufts master

